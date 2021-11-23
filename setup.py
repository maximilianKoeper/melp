import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setuptools.setup(
    install_requires=required,
    name="melp",
    version="0.0.1",
    author="Maximilian Koeper / Erik Steinkamp",
    author_email="maximilian.koeper@kip.uni-heidelberg.de",
    description="Mu3eHelper for tile analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maximilianKoeper/melp",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "melp"},
    packages=setuptools.find_packages(where="melp"),
    python_requires=">=3.6",
)