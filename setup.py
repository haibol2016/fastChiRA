#!/usr/bin/env python
"""Setup script for fastChiRA. Prefer installing via: pip install .

Primary metadata is also declared in pyproject.toml for modern tooling.
"""

import os
import re

from setuptools import find_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))


def read_version():
    """Read __version__ from chira.utilities (single source of truth)."""
    version_file = os.path.join(HERE, "src", "chira", "utilities.py")
    with open(version_file, encoding="utf-8") as f:
        for line in f:
            match = re.search(
                r"^__version__\s*=\s*['\"]([\d]+\.[\d]+\.[\d]+)['\"]", line
            )
            if match:
                return match.group(1)
    return "0.0.0"


def read_readme():
    readme = os.path.join(HERE, "README.md")
    if os.path.exists(readme):
        with open(readme, encoding="utf-8") as f:
            return f.read()
    return ""


setup(
    name="chira",
    version=read_version(),
    description="fastChiRA - Chimeric Read Analyzer for RNA-RNA interactome data",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    author="fastChiRA Contributors",
    url="https://github.com/haibol2016/chira",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={
        "chira": ["share/*.R", "share/*.tmpl"],
    },
    include_package_data=True,
    python_requires=">=3.6",
    install_requires=[
        "biopython",
        "bcbio-gff",
        "pysam",
        "mpire>=2.4.0",
    ],
    extras_require={
        "optional": [
            "psutil",
            "requests",
            "pyliftover",
        ],
    },
    entry_points={
        "console_scripts": [
            "chira_collapse.py=chira.collapse:main",
            "chira_map.py=chira.map:main",
            "chira_merge.py=chira.merge:main",
            "chira_quantify.py=chira.quantify:main",
            "chira_extract.py=chira.extract:main",
            "concatenate_gtf.py=chira.utils.concatenate_gtf:main",
            "download_ensembl.py=chira.utils.download_ensembl:main",
            "download_mirbase_gff3.py=chira.utils.download_mirbase_gff3:main",
            "download_mirbase_mature.py=chira.utils.download_mirbase_mature:main",
            "extract_transcripts_from_genome.py=chira.utils.extract_transcripts_from_genome:main",
            "remove_mirna_hairpin_from_gtf.py=chira.utils.remove_mirna_hairpin_from_gtf:main",
            "process_map_chunk_batchtools.py=chira.batchtools.process_map_chunk:main",
            "process_intarna_chunk_batchtools.py=chira.batchtools.process_intarna_chunk:main",
            "merge_intarna_into_chimeras.py=chira.batchtools.merge_intarna:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
