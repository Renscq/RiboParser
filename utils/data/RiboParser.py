#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : riboparser.py


import pkg_resources


class RiboParserInfo:
    VERSION = "0.1.14"
    UPDATE_DATE = "2024-06-10"
    CITATION = (
        '''
        Shuchao Ren, Yinan Li, Zhipeng Zhou. 
        RiboParser/RiboShiny: An integrated platform for comprehensive analysis and visualization of ribo-seq data. 
        Journal of Genetics and Genomics (2025) 
        doi:10.1016/j.jgg.2025.04.010.
        '''
    )
    REQUIRED_PACKAGES = ["pandas", "polars", "numpy", "matplotlib-venn", "seqlogo", 
                         "matplotlib", "seaborn", "biopython", 
                         "scipy", "scikit-learn", "statsmodels", 
                         "pysam", "joblib"]

    @classmethod
    def show_version(cls):
        print(f"RiboParser version: {cls.VERSION}")
        print(f"Last update: {cls.UPDATE_DATE}")

    @classmethod
    def show_citation(cls):
        print("Please cite:")
        print(cls.CITATION)

    @classmethod
    def check_dependencies(cls):
        missing = []
        for pkg in cls.REQUIRED_PACKAGES:
            try:
                pkg_resources.get_distribution(pkg)
            except pkg_resources.DistributionNotFound:
                missing.append(pkg)
        if missing:
            print(f"Missing dependencies: {', '.join(missing)}")
            return False
        else:
            print(cls.REQUIRED_PACKAGES)
        print("All required dependencies are installed.")
        return True

    @classmethod
    def check_package_modules(cls):

        from pathlib import Path

        project_path = Path(__file__).resolve()
        # check the root directory（directory contains include pyproject.toml / README.md / .git）
        for _ in range(8):
            if any((project_path / name).exists() for name in ("pyproject.toml", "README.md", ".git")):
                break
            if project_path.parent == project_path:
                break
            project_path = project_path.parent
        root = project_path

        utils_dir = root / "utils"
        scripts_dir = root / "scripts"

        @staticmethod
        def module_name_from_path(p: Path):
            try:
                rel = p.relative_to(root)
            except Exception:
                rel = p
            return ".".join(rel.with_suffix("").parts)

        rpf = []
        serp = []
        classes = []
        others = []

        if utils_dir.exists():
            for now_path in utils_dir.iterdir():
                if now_path.is_file() and now_path.suffix == ".py" and not now_path.name.startswith("_"):
                    mod = module_name_from_path(now_path)
                    name = now_path.stem
                    if name.startswith("rpf_") or name.startswith("rna_"):
                        rpf.append(mod)
                    elif name.startswith("serp_"):
                        serp.append(mod)
                    else:
                        others.append(mod)
                elif now_path.is_dir():
                    for sub in now_path.rglob("*.py"):
                        if sub.name.startswith("_") or sub.name == "__init__.py":
                            continue
                        mod = module_name_from_path(sub)
                        if sub.stem.startswith("rpf_"):
                            rpf.append(mod)
                        elif sub.stem.startswith("serp_"):
                            serp.append(mod)
                        else:
                            classes.append(mod)

        if scripts_dir.exists():
            for now_path in scripts_dir.rglob("*.py"):
                if now_path.name.startswith("_") or now_path.name == "__init__.py":
                    continue
                mod = module_name_from_path(now_path)
                if now_path.stem.startswith("rpf_"):
                    rpf.append(mod)
                elif now_path.stem.startswith("serp_"):
                    serp.append(mod)
                else:
                    others.append(mod)

        # 去重并排序
        rpf = sorted(set(rpf))
        serp = sorted(set(serp))
        classes = sorted(set(classes))
        others = sorted(set(others))

        @staticmethod
        def try_import(module_name: str) -> bool:
            try:
                import importlib

                importlib.import_module(module_name)
                return True
            except Exception:
                return False

        print("RPF modules:")
        if rpf:
            for now_module in rpf:
                status = "[import OK]" if try_import(now_module) else "[import FAILED]"
                print(f" - {now_module} {status}")
        else:
            print(" - (not found)")

        print("SERP modules:")
        if serp:
            for now_module in serp:
                status = "[import OK]" if try_import(now_module) else "[import FAILED]"
                print(f" - {now_module} {status}")
        else:
            print(" - (not found)")

        print("Classes:")
        if classes:
            for now_module in classes:
                status = "[import OK]" if try_import(now_module) else "[import FAILED]"
                print(f" - {now_module} {status}")
        else:
            print(" - (not found)")

        print("Other scripts:")
        if others:
            for now_module in others:
                status = "[import OK]" if try_import(now_module) else "[import FAILED]"
                print(f" - {now_module} {status}")
        else:
            print(" - (not found)")
