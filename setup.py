import os
import setuptools
from distutils.core import setup


def get_motevo_files():
    path = os.path.join("motevowrapper", "tests", "data")
    matches = os.listdir(path)
    matches = [os.path.join(path, m) for m in matches]
    return matches


def get_version():
    with open("motevowrapper/_version.py", "r") as f:
        for line in f:
            if line.startswith("version"):
                return line.split("=")[1].strip()[1:-1]


version = get_version()

with open("README.md", "r", encoding="utf-8") as f:
    README = f.read()

setup(
    name="motevowrapper",
    packages=["motevowrapper"],
    version=version,
    license="MIT",
    description="Simple Python MotEvo wrapper.",
    author="Đorđe Relić",
    author_email="dorde.relic@protonmail.com",
    url="https://github.com/brlauuu/motevowrapper",
    download_url=f"https://github.com/brlauuu/motevowrapper/archive/v{version}.tar.gz",
    keywords=["MotEvo", "wrapper", "binding", "sites", "tfbs-discovery"],
    install_requires=["pandas",],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    long_description_content_type="text/markdown",
    long_description=README,
    include_package_data=True,
    package_data={"motevowrapper": get_motevo_files(),},
)
