from setuptools import setup

setup(
    name='riboParser',
    version='1.0',
    packages=['scripts', 'scripts.foo'],
    install_requires=['pysam', 'matplotlib', 'pandas', 'numpy', 'biopython'],

    entry_points={

        'console_scripts': [
            # name determined the name of cmd line direct call
            'make_ensb_ref=scripts.make_ensb_ref:main',
            'detect_offset=scripts.detect_offset:main',
            'ribo_parser=scripts.ribo_parser:main',
            'meta_plot=scripts.meta_plot:main',
            'rpf_plot=scripts.rpf_plot:main',
            'cdf_plot=scripts.cdf_plot:main',
            'cdt_calc=scripts.cdt_calc:main'
        ],
    },

    url='',
    license='MIT',
    author='Ren Shuchao',
    author_email='rensc0718@163.com',
    description='A toolkit for ribosome sequence parsing.'
)
