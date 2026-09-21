import os
import subprocess
import argparse
import sys
from pathlib import Path


class LinkerDesignPipeline:
    def __init__(self, args):
        self.args = args

        # --- 源码/工具目录 ---
        self.root_tool_dir = "/root/autodl-fs"
        self.rfdiffusion_src = os.path.join(self.root_tool_dir, "RFdiffusion")
        self.mpnn_src = os.path.join(self.root_tool_dir, "ProteinMPNN")
        self.analysis_script = "/root/autodl-fs/af2_results/analyze_linker.py"

        # --- 统一输出目录结构 ---
        self.base_dir = Path(args.output_dir).resolve()
        self.rfd_out_dir = self.base_dir / "01_rfdiffusion"
        self.mpnn_out_dir = self.base_dir / "02_mpnn"
        self.af2_out_dir = self.base_dir / "03_colabfold"
        self.analysis_out_dir = self.base_dir / "04_analysis"

        # 创建目录结构
        for d in [self.rfd_out_dir, self.mpnn_out_dir, self.af2_out_dir, self.analysis_out_dir]:
            d.mkdir(parents=True, exist_ok=True)

        print(f"[INIT] Output directory set to: {self.base_dir}")

    def run_command(self, cmd, cwd=None):
        """执行Shell命令并打印日志"""
        print(f"\n[INFO] Running: {cmd}")
        try:
            subprocess.check_call(cmd, shell=True, cwd=cwd, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Command failed with error code {e.returncode}")
            sys.exit(1)

    def step_1_rfdiffusion(self):
        """Step 1: RFdiffusion (Backbone Generation)"""
        print("### Step 1: RFdiffusion ###")

        contig_str = (
            f"'{self.args.chain}{self.args.phage_start}-{self.args.phage_end}/"
            f"{self.args.linker_min}-{self.args.linker_max}/"
            f"{self.args.chain}{self.args.pep_start}-{self.args.pep_end}'"
        )

        output_prefix = self.rfd_out_dir / "design"

        cmd = (
            F"run_inference.py "
            f"inference.input_pdb={self.args.pdb_path} "
            f"inference.output_prefix={output_prefix} "
            f"inference.num_designs={self.args.num_designs} "
            f"contigmap.contigs=[{contig_str}]"
        )

        self.run_command(cmd, cwd=self.rfdiffusion_src)
        return self.rfd_out_dir

    def step_2_protein_mpnn(self, rfdiffusion_out_dir):
        """Step 2: ProteinMPNN (Sequence Design)"""
        print("### Step 2: ProteinMPNN ###")

        phage_len = self.args.phage_end - self.args.phage_start + 1
        pep_len = self.args.pep_end - self.args.pep_start + 1

        fixed_pos_file = self.mpnn_out_dir / "fixed_positions.jsonl"
        parsed_batch_file = self.mpnn_out_dir / "parsed_batch.jsonl"
        mpnn_seqs_folder = self.mpnn_out_dir / "final_seqs"

        capsid_len_param = self.args.capsid_len_override if self.args.capsid_len_override else phage_len

        cmd_fix = (
            f"python make_fixed_positions_batch.py "
            f"--input_folder '{rfdiffusion_out_dir}' "
            f"--capsid_len {capsid_len_param} "
            f"--avitag_len {pep_len} "
            f"--output_path '{fixed_pos_file}'"
        )
        self.run_command(cmd_fix, cwd=self.mpnn_src)

        cmd_parse = (
            f"python helper_scripts/parse_multiple_chains.py "
            f"--input_path='{rfdiffusion_out_dir}' "
            f"--output_path='{parsed_batch_file}'"
        )
        self.run_command(cmd_parse, cwd=self.mpnn_src)

        cmd_run = (
            f"python protein_mpnn_run.py "
            f"--jsonl_path '{parsed_batch_file}' "
            f"--fixed_positions '{fixed_pos_file}' "
            f"--out_folder '{mpnn_seqs_folder}' "
            f"--num_seq_per_target {self.args.seqs_per_design} "
            f"--sampling_temp 0.2 "
            f"--batch_size 1"
        )
        self.run_command(cmd_run, cwd=self.mpnn_src)

        return mpnn_seqs_folder / "seqs"

    def step_3_colabfold(self, mpnn_seqs_dir):
        """Step 3: ColabFold (Validation) - Loop Mode with Subdirectories"""
        print("### Step 3: ColabFold ###")

        input_dir = Path(mpnn_seqs_dir)
        fasta_files = sorted(list(input_dir.glob("*.fa")))

        if not fasta_files:
            raise FileNotFoundError(f"No fasta files found in {input_dir}")

        total_files = len(fasta_files)
        print(f"[INFO] Found {total_files} FASTA files. Running ColabFold sequentially...")

        for index, fasta_file in enumerate(fasta_files, 1):
            print(f"\n--- [Task {index}/{total_files}] Processing: {fasta_file.name} ---")

            # --- 修改开始：创建独立的子文件夹 ---
            # 使用文件名（如 design_0）作为子目录名
            design_name = fasta_file.stem  # 获取文件名不带后缀，例如 "design_0"
            current_output_dir = self.af2_out_dir / design_name
            current_output_dir.mkdir(parents=True, exist_ok=True)
            # --- 修改结束 ---

            # 将输出目录指向子文件夹
            cmd = (
                f"colabfold_batch '{fasta_file}' '{current_output_dir}' "
                f"--num-recycle 3 --use-gpu-relax"
            )

            self.run_command(cmd, cwd=self.base_dir)

        return self.af2_out_dir

    def step_4_analysis(self, af2_result_dir):
        """Step 4: Analysis (Auto-Length Detection) - Supports Subdirectories"""
        print("### Step 4: Analysis ###")

        phage_len = self.args.phage_end - self.args.phage_start + 1
        pep_len = self.args.pep_end - self.args.pep_start + 1

        print(f"[INFO] Auto-Detect Params: Phage={phage_len}, Peptide={pep_len}")

        input_root = Path(af2_result_dir)

        # 检查是否存在子目录（判断是否使用了修改版的 Step 3）
        subdirs = sorted([d for d in input_root.iterdir() if d.is_dir()])

        if subdirs:
            print(f"[INFO] Detected {len(subdirs)} subdirectories. Running analysis for each design...")

            for sub in subdirs:
                design_name = sub.name  # 例如 "design_0"

                # 为每个 design 创建独立的分析输出目录，防止结果覆盖
                current_analysis_out = self.analysis_out_dir / design_name
                current_analysis_out.mkdir(parents=True, exist_ok=True)

                print(f"   -> Analyzing: {design_name}")

                cmd = (
                    f"python {self.analysis_script} "
                    f"--input '{sub}' "
                    f"--phage_len {phage_len} "
                    f"--pep_len {pep_len} "
                    f"--output '{current_analysis_out}'"
                )

                # 这里捕获错误但不退出，防止一个分析失败导致整体中断
                try:
                    subprocess.check_call(cmd, shell=True, executable='/bin/bash', cwd=self.root_tool_dir)
                except subprocess.CalledProcessError:
                    print(f"[WARNING] Analysis failed for {design_name}, skipping...")

        else:
            # 如果没有子目录（即 Step 3 还是扁平输出），则保持原有逻辑
            print("[INFO] No subdirectories found. Running analysis on root folder...")
            cmd = (
                f"python {self.analysis_script} "
                f"--input '{af2_result_dir}' "
                f"--phage_len {phage_len} "
                f"--pep_len {pep_len} "
                f"--output '{self.analysis_out_dir}'"
            )
            self.run_command(cmd, cwd=self.root_tool_dir)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Automated Linker Design Pipeline")

    parser.add_argument("--pdb_path", required=True, help="Absolute path to input PDB file")
    parser.add_argument("--output_dir", default="./unified_output", help="Root directory for all pipeline outputs")
    parser.add_argument("--chain", default="A", help="Chain ID in the PDB")
    parser.add_argument("--phage_start", type=int, required=True, help="Start residue of Phage")
    parser.add_argument("--phage_end", type=int, required=True, help="End residue of Phage")
    parser.add_argument("--linker_min", type=int, default=10, help="Min length of linker")
    parser.add_argument("--linker_max", type=int, default=15, help="Max length of linker")
    parser.add_argument("--pep_start", type=int, required=True, help="Start residue of Peptide")
    parser.add_argument("--pep_end", type=int, required=True, help="End residue of Peptide")
    parser.add_argument("--num_designs", type=int, default=10, help="Number of RFdiffusion designs")
    parser.add_argument("--seqs_per_design", type=int, default=2, help="Sequences per structure in MPNN")
    parser.add_argument("--capsid_len_override", type=int, help="Override capsid_len for MPNN if needed")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    pipeline = LinkerDesignPipeline(args)

    rfd_out = pipeline.step_1_rfdiffusion()
    mpnn_seqs = pipeline.step_2_protein_mpnn(rfd_out)
    af2_out = pipeline.step_3_colabfold(mpnn_seqs)
    pipeline.step_4_analysis(af2_out)

    print(f"\n[SUCCESS] Pipeline completed. All results in: {args.output_dir}")
