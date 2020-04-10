from setuptools import setup

requirements = [
    'scikit-allel',
    'numpy',
    'pandas',
	'dill',
]

setup(
    name='pixy',
    version='0.94.04',
    packages=['pixy'],
    entry_points={
        'console_scripts': [
            'pixy=pixy.__main__:main'
        ]
    },
    url='https://github.com/ksamuk/pixy',
    license='MIT',
    author='Katharine Korunes, Kieran Samuk',
    author_email='kkorunes@gmail.com, ksamuk@gmail.com',
    description="pixy",
    install_requires=requirements,
    keywords='pixy',
    classifiers=[
        'Programming Language :: Python :: 3.6'
    ]
)