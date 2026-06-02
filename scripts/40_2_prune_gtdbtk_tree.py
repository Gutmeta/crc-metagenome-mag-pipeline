#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path

def strip_newick_comments(newick_text: str) -> str:
    """Remove Newick comments like [&&NHX:...] or [....] which may break parsers."""
    return re.sub(r"\[[^\]]*\]", "", newick_text)

def norm_label(x: str, strip_ext: bool) -> str:
    if not strip_ext:
        return x
    return re.sub(r"\.(fa|fna|fasta)$", "", x)

def is_leaf_compat(node) -> bool:
    """
    ete versions differ:
      - ete3: node.is_leaf() method
      - some ete4 builds: node.is_leaf bool attribute
    """
    v = getattr(node, "is_leaf", None)
    if callable(v):
        return v()
    if isinstance(v, bool):
        return v
    return len(getattr(node, "children", [])) == 0

def leaf_names_compat(t):
    """
    Return leaf names as a Python list across different ete4 builds.
    """
    # preferred
    if hasattr(t, "leaf_names"):
        ln = t.leaf_names()
        return list(ln)  # works if ln is list or generator
    # fallbacks
    if hasattr(t, "get_leaf_names"):
        return list(t.get_leaf_names())
    if hasattr(t, "iter_leaf_names"):
        return list(t.iter_leaf_names())
    # last resort: traverse
    return [n.name for n in t.traverse() if is_leaf_compat(n)]

def main():
    ap = argparse.ArgumentParser(
        description="Prune a Newick tree to keep only selected tips (e.g., C1A__/C1B__), "
                    "and optionally remove internal node labels (GTDB decorated labels like f__/o__)."
    )
    ap.add_argument("-i", "--input_tree",
                    default="./gtdbtk_denovo/gtdbtk.bac120.decorated.tree",
                    help="Input tree file (Newick).")
    ap.add_argument("-o", "--output_tree",
                    default="./gtdbtk_denovo/tree.C1AB.nwk",
                    help="Output pruned tree file (Newick).")
    ap.add_argument("--pattern",
                    default=r"^C1[AB]__",
                    help=r"Regex pattern to match tip labels to KEEP. Default: ^C1[AB]__")
    ap.add_argument("--keep_list",
                    default=None,
                    help="Optional file with tip IDs to keep (one per line). Overrides --pattern.")
    ap.add_argument("--strip_ext",
                    action="store_true",
                    help="If set, match IDs after stripping .fa/.fna/.fasta extensions.")
    ap.add_argument("--remove_internal_labels",
                    action="store_true",
                    default=True,
                    help="If set, remove internal node labels after pruning (recommended for GTDB decorated trees).")
    ap.add_argument("--parser",
                    default="name",
                    help="ETE newick parser. For GTDB decorated trees, use 'name' (default).")
    ap.add_argument("--check_only",
                    action="store_true",
                    help="Only report how many tips match, do not write output.")
    args = ap.parse_args()

    in_path = Path(args.input_tree)
    out_path = Path(args.output_tree)

    # Read and sanitize Newick
    newick = in_path.read_text(encoding="utf-8", errors="ignore").strip()
    if not newick.endswith(";"):
        newick += ";"
    newick = strip_newick_comments(newick)

    # Import ete4
    try:
        from ete4 import Tree
    except ImportError:
        print(
            "ERROR: ete4 is not installed.\n"
            "Install with one of:\n"
            "  pip install ete4\n"
            "  conda install -c etetoolkit ete4\n",
            file=sys.stderr
        )
        sys.exit(1)

    # Parse tree (GTDB decorated trees: parser='name' is safest)
    try:
        t = Tree(newick, parser=args.parser)
    except Exception as e:
        print(f"ERROR: Failed to parse tree with parser={args.parser}: {e}", file=sys.stderr)
        sys.exit(2)

    # Get leaf names as list
    leaf_names = leaf_names_compat(t)

    # Determine tips to keep
    keep_names = []
    if args.keep_list:
        wanted = []
        with open(args.keep_list, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    wanted.append(line)
        wanted_set = set(norm_label(x, args.strip_ext) for x in wanted)
        for name in leaf_names:
            if norm_label(name, args.strip_ext) in wanted_set:
                keep_names.append(name)
    else:
        pat = re.compile(args.pattern)
        for name in leaf_names:
            if pat.search(norm_label(name, args.strip_ext)):
                keep_names.append(name)

    if not keep_names:
        print("ERROR: No tips matched. Check --pattern/--keep_list and whether tip labels match.", file=sys.stderr)
        print("Sample tip labels from tree:", file=sys.stderr)
        for s in leaf_names[:20]:
            print("  ", s, file=sys.stderr)
        sys.exit(3)

    print(f"[INFO] Total tips in input tree: {len(leaf_names)}")
    print(f"[INFO] Tips matched for pruning: {len(keep_names)}")

    if args.check_only:
        nonmatch = [x for x in leaf_names if x not in set(keep_names)]
        print(f"[INFO] Non-matching tips (first 10): {nonmatch[:10]}")
        return

    # Prune
    t.prune(keep_names, preserve_branch_length=True)

    # Remove internal node labels if requested (keeps leaf names intact)
    if args.remove_internal_labels:
        for n in t.traverse():
            if not is_leaf_compat(n):
                n.name = ""

    # Write output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        t.write(outfile=str(out_path), parser="name")
    except Exception as e:
        print(f"ERROR: Failed to write pruned tree: {e}", file=sys.stderr)
        sys.exit(4)

    # Post-check
    try:
        t2 = Tree(out_path.read_text(encoding="utf-8", errors="ignore"), parser="name")
        tips2 = leaf_names_compat(t2)
        bad = [x for x in tips2 if not (x.startswith("C1A__") or x.startswith("C1B__"))]
        print(f"[INFO] Output tips: {len(tips2)}; non-C1 tips: {len(bad)}")
        if bad:
            print("[WARN] Examples of non-C1 tips in output:", bad[:10])
    except Exception:
        pass

    print(f"[INFO] Wrote pruned tree: {out_path}")

if __name__ == "__main__":
    main()