from setuptools import setup

requirements = [
    'scikit-allel',
    'tqdm'
]

setup(
    name='pixy',
    version='1.0',
    description="pixy test",
    #author="kieran",
    #author_email='ksamuk',
    url='https://github.com/ksamuk/pixy',
    packages=['pixy'],
    entry_points={
        'console_scripts': [
            'pixy=pixy.py'
        ]
    },
	scripts='src/pixy.py',
    install_requires=requirements,
    keywords='pixy',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'
    ]
)
