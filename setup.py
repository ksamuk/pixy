from setuptools import setup

requirements = ["scikit-allel", "pandas", "numpy", "multiprocess", "scipy", "numcodecs"]

setup(
    name="pixy",
    version="2.0.0.beta1",
    packages=["pixy"],
    entry_points={"console_scripts": ["pixy=pixy.__main__:main"]},
    url="https://github.com/ksamuk/pixy",
    license="MIT",
    author="Kieran Samuk",
    author_email="ksamuk@gmail.com",
    description="pixy",
    install_requires=requirements,
    keywords="pixy",
    classifiers=["Programming Language :: Python :: 3.9"],
)
