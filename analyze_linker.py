import os
import glob
import csv
import argparse
import warnings
import numpy as np

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio import SeqUtils

    warnings.simplefilter("ignore", PDBConstructionWarning)
except ImportError:
    raise SystemExit("错误: 请先安装 biopython (pip install biopython)")

def _is_protein_residue(res):
    # 只保留标准蛋白残基，去掉水/配体等
    return res.id[0] == " "

def extract_plddt_from_residue(res):
    """
    AlphaFold PDB 通常把 pLDDT 放在每个原子的 B-factor。
    优先取 CA 的 B-factor；没有 CA 就取该残基所有原子 B-factor 的平均。
    """
    try:
        if "CA" in res:
            return float(res["CA"].bfactor)
    except Exception:
        pass

    bfs = []
    for atom in res.get_atoms():
        try:
            bfs.append(float(atom.bfactor))
        except Exception:
            continue
    return float(np.mean(bfs)) if bfs else None

def get_structural_metrics_from_pdb(pdb_path, start_idx, end_idx, chain_id=None):
    """
    Returns:
      dist (float),
      avg_sasa (float),
      seq (str),
      avg_plddt (float or None)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    model = structure[0]

    # 选链：默认第一条链；也可用 --chain 指定
    chain = None
    if chain_id:
        if chain_id in model:
            chain = model[chain_id]
        else:
            return None, None, "", None
    else:
        chains = list(model.get_chains())
        if not chains:
            return None, None, "", None
        chain = chains[0]

    residues_all = [r for r in chain.get_residues() if _is_protein_residue(r)]
    n_total = len(residues_all)
    if n_total == 0:
        return None, None, "", None

    safe_start = max(0, min(start_idx, n_total - 1))
    safe_end = max(0, min(end_idx, n_total))
    if safe_start >= safe_end:
        return 0.0, 0.0, "", None

    # 1) Linker 序列
    seq = ""
    for i in range(safe_start, safe_end):
        res = residues_all[i]
        try:
            seq += SeqUtils.seq1(res.get_resname())
        except Exception:
            seq += "X"

    # 2) End-to-End Distance (CA-CA)
    dist = 0.0
    try:
        n_res = residues_all[safe_start]
        c_res = residues_all[safe_end - 1]
        if "CA" in n_res and "CA" in c_res:
            diff = n_res["CA"].coord - c_res["CA"].coord
            dist = float(np.sqrt(np.sum(diff ** 2)))
    except Exception:
        dist = 0.0

    # 3) SASA（残基层级）
    avg_sasa = 0.0
    try:
        sr = ShrakeRupley()
        sr.compute(model, level="R")
        sasa_vals = []
        for i in range(safe_start, safe_end):
            res = residues_all[i]
            if hasattr(res, "sasa"):
                sasa_vals.append(float(res.sasa))
        avg_sasa = float(np.mean(sasa_vals)) if sasa_vals else 0.0
    except Exception:
        avg_sasa = 0.0

    # 4) 平均 pLDDT（从 B-factor 读）
    plddts = []
    for i in range(safe_start, safe_end):
        v = extract_plddt_from_residue(residues_all[i])
        if v is not None:
            plddts.append(v)
    avg_plddt = float(np.mean(plddts)) if plddts else None

    return dist, avg_sasa, seq, avg_plddt

def evaluate_linker_tags(avg_plddt, dist, sasa):
    """
    沿用你原来的判定阈值（去掉 PAE 相关的 Rigid/Flexible）。
    """
    tags = []

    if avg_plddt is not None:
        if avg_plddt > 80:
            tags.append("Stable")
        elif avg_plddt < 50:
            tags.append("Disordered")
    else:
        tags.append("pLDDT_NA")

    if sasa > 40:
        tags.append("Exposed")
    elif sasa < 20:
        tags.append("Buried")

    if dist < 5.0:
        tags.append("Collapsed")
    elif dist > 15.0:
        tags.append("Extended")

    return " ".join(tags) if tags else "Unknown"

def batch_analyze_linker_pdb_only(input_folder, phage_len, pep_len, output_dir="./linker_results", chain_id=None):
    if not os.path.exists(input_folder):
        print(f"错误: 输入文件夹不存在 -> {input_folder}")
        return

    pdb_files = []
    # 常见：*_unrelaxed_*.pdb 或任意 .pdb
    pdb_files += glob.glob(os.path.join(input_folder, "*unrelaxed*.pdb"))
    pdb_files += glob.glob(os.path.join(input_folder, "*.pdb"))
    pdb_files = sorted(list(dict.fromkeys(pdb_files)))  # 去重保持顺序

    if not pdb_files:
        print(f"错误: 在 {input_folder} 没有找到 .pdb 文件。")
        return

    os.makedirs(output_dir, exist_ok=True)
    results = []

    print(f"正在分析 {len(pdb_files)} 个 PDB...")
    print(f"固定区域设置 -> Phage前段长度: {phage_len}, Peptide后段长度: {pep_len}")
    if chain_id:
        print(f"指定链 -> {chain_id}")
    print(f"结果将保存至 -> {os.path.join(output_dir, 'summary_results.csv')}\n")

    for pdb_path in pdb_files:
        filename = os.path.basename(pdb_path)
        sample_name = os.path.splitext(filename)[0]

        # 需要先知道总长度：这里直接解析一次拿 residues 数
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("model", pdb_path)
            model = structure[0]
            chain = model[chain_id] if chain_id else list(model.get_chains())[0]
            residues_all = [r for r in chain.get_residues() if _is_protein_residue(r)]
            total_len = len(residues_all)
        except Exception:
            total_len = 0

        if total_len <= 0:
            results.append({
                "Sample Name": sample_name,
                "Sequence": "",
                "pLDDT": "N/A",
                "PAE": "N/A",
                "Dist_A": "N/A",
                "SASA": "N/A",
                "Status": "⚠️ Err",
                "Evaluation": "ParseFailed"
            })
            continue

        linker_start_idx = phage_len
        linker_end_idx = total_len - pep_len

        if linker_end_idx <= linker_start_idx:
            print(f"[Warn] Sample {sample_name} too short (total {total_len}). Skipping.")
            continue

        dist, sasa, seq, avg_plddt = get_structural_metrics_from_pdb(
            pdb_path, linker_start_idx, linker_end_idx, chain_id=chain_id
        )

        if dist is None:
            results.append({
                "Sample Name": sample_name,
                "Sequence": "",
                "pLDDT": "N/A",
                "PAE": "N/A",
                "Dist_A": "N/A",
                "SASA": "N/A",
                "Status": "⚠️ Err",
                "Evaluation": "ParseFailed"
            })
            continue

        evaluation = evaluate_linker_tags(avg_plddt, dist, sasa)

        results.append({
            "Sample Name": sample_name,
            "Sequence": seq,
            "pLDDT": round(avg_plddt, 2) if avg_plddt is not None else "N/A",
            "PAE": "N/A",
            "Dist_A": round(dist, 2),
            "SASA": round(sasa, 2),
            "Status": "✅ OK",
            "Evaluation": evaluation
        })

    if not results:
        print("没有产生有效分析结果。")
        return

    # 排序：优先 Stable/Exposed 之类不太好量化；这里按 pLDDT(降序) + SASA(降序) + Dist(降序) 排
    def _plddt_key(x):
        return float(x["pLDDT"]) if isinstance(x["pLDDT"], (int, float)) else -1.0

    results.sort(key=lambda x: (_plddt_key(x), float(x["SASA"]) if x["SASA"] != "N/A" else -1.0, float(x["Dist_A"]) if x["Dist_A"] != "N/A" else -1.0), reverse=True)

    csv_file = os.path.join(output_dir, "summary_results.csv")
    csv_columns = ["Sample Name", "Sequence", "pLDDT", "PAE", "Dist_A", "SASA", "Status", "Evaluation"]

    with open(csv_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"[SUCCESS] 结果已保存至: {csv_file}\n")

    print("-" * 130)
    print(f"{'Sample Name':<28} | {'Sequence':<18} | {'pLDDT':<6} | {'PAE':<6} | {'Dist(Å)':<7} | {'SASA':<6} | {'Eval'}")
    print("-" * 130)
    for res in results:
        name = res["Sample Name"][:26]
        seq = res["Sequence"]
        seq_disp = (seq[:15] + "...") if len(seq) > 15 else seq
        print(f"{name:<28} | {seq_disp:<18} | {res['pLDDT']:<6} | {res['PAE']:<6} | {res['Dist_A']:<7} | {res['SASA']:<6} | {res['Evaluation']}")
    print("-" * 130)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate linker effect using PDB-only metrics (Dist/SASA/pLDDT-from-Bfactor).")
    parser.add_argument("--input", "-i", type=str, required=True, help="Folder containing PDB files")
    parser.add_argument("--phage_len", type=int, required=True, help="Length of fixed Phage region (N-terminus)")
    parser.add_argument("--pep_len", type=int, required=True, help="Length of fixed Peptide region (C-terminus)")
    parser.add_argument("--chain", type=str, default=None, help="Chain ID to use (optional). Default: first chain")
    parser.add_argument("--output", "-o", type=str, default="./linker_results", help="Output directory for CSV")

    args = parser.parse_args()
    batch_analyze_linker_pdb_only(args.input, args.phage_len, args.pep_len, args.output, chain_id=args.chain)
