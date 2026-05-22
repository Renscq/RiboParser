#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : riboparser.py


from importlib.metadata import version, PackageNotFoundError


class RiboParserInfo:
    try:
        version = version("RiboParser")
    except PackageNotFoundError:
        version = "unknown"
    
    update_date = "2026-05-21"
    citation = (
        '''
        Shuchao Ren, Yinan Li, Zhipeng Zhou. 
        RiboParser/RiboShiny: An integrated platform for comprehensive analysis and visualization of ribo-seq data. 
        Journal of Genetics and Genomics (2025) 
        doi:10.1016/j.jgg.2025.04.010.
        '''
    )
    required_packages = ["pandas", "polars", "numpy", "matplotlib-venn", "seqlogo", 
                         "matplotlib", "seaborn", "biopython", 
                         "scipy", "scikit-learn", "statsmodels", 
                         "pysam", "joblib"]

    @classmethod
    def show_version(cls):
        print(f"RiboParser version: {cls.version}")
        print(f"Last update: {cls.update_date}")

    @classmethod
    def show_citation(cls):
        print("Please cite:")
        print(cls.citation)

    @classmethod
    def check_dependencies(cls):
        missing = []
        for pkg in cls.required_packages:
            try:
                version(pkg)
            except PackageNotFoundError:
                missing.append(pkg)
        if missing:
            print(f"Missing dependencies: {', '.join(missing)}")
            return False
        else:
            print(cls.required_packages)
        print("All required dependencies are installed.")
        return True

    @classmethod
    def check_package_modules(cls, module_type: str = "all"):
        from pathlib import Path
        import sys
        import importlib

        script_path = Path(__file__).resolve()

        # Find project root
        root = script_path.parent
        for _ in range(10):
            if any((root / name).exists() for name in ("pyproject.toml", "README.md", ".git", "utils", "scripts")):
                break
            if root.parent == root:
                break
            root = root.parent

        # Make local modules importable
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))

        utils_dir = root / "utils"
        scripts_dir = root / "scripts"

        modules = {
            "ribo": [],
            "serp": [],
            "smorf": [],
            "scripts": []
        }

        def module_name_from_path(p: Path) -> str:
            rel = p.relative_to(root)
            return ".".join(rel.with_suffix("").parts)

        def add_module(p: Path):
            if p.name.startswith("_") or p.name == "__init__.py":
                return

            mod = module_name_from_path(p)
            parts = p.relative_to(root).parts
            stem = p.stem

            if "smorf" in parts or stem.startswith("smorf_"):
                modules["smorf"].append(mod)
            elif "serp" in parts or stem.startswith("serp_"):
                modules["serp"].append(mod)
            elif "ribo" in parts or stem.startswith(("rpf_", "rna_")):
                modules["ribo"].append(mod)
            elif "scripts" in parts:
                modules["scripts"].append(mod)

        if utils_dir.exists():
            for p in utils_dir.rglob("*.py"):
                add_module(p)

        if scripts_dir.exists():
            for p in scripts_dir.rglob("*.py"):
                add_module(p)

        for key in modules:
            modules[key] = sorted(set(modules[key]))

        def try_import(module_name: str) -> bool:
            try:
                importlib.import_module(module_name)
                return True
            except Exception as e:
                return False

        show_keys = modules.keys() if module_type == "all" else [module_type]

        for key in show_keys:
            print(f"{key} modules:")
            if modules.get(key):
                for mod in modules[key]:
                    status = "[import OK]" if try_import(mod) else "[import FAILED]"
                    print(f" - {mod} {status}")
            else:
                print(" - (not found)")
