from setuptools import setup

requirements = [
    'scikit-allel',
    'tqdm',
    'numpy',
    'pandas',
]

setup(
    name='pixy',
    version='0.1',
    packages=[''],
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'pixy=pixy.py'
        ]
    },
    scripts="pixy.py",
    url='https://github.com/ksamuk/pixy',
    license='MIT',
    author='Katharine Korunes, Kieran Samuk',
    author_email='kkorunes@gmail.com, ksamuk@gmail.com',
    description="pixy",
    install_requires=requirements,
    keywords='pixy',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'
    ]
)
