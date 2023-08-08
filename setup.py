from setuptools import setup, find_packages

with open('requirements.txt') as f:
	requirements = f.readlines()

long_description = 'Create interactive and customizable gene transcript visualizations'

setup(
		name ='transcriptionary',
		version ='0.1.1',
		author ='Suchita Lulla',
		author_email ='suchitalulla@gmail.com',
		url ='https://github.com/quinlan-lab/transcriptionary',
		description ='transcript plots',
		long_description = long_description,
		long_description_content_type ="text/markdown",
		license ='MIT',
		packages = find_packages(),
		entry_points ={
			'console_scripts': [
                'transcriptionary = transcriptionary.transcriptionary:main'
			]
		},
		classifiers =(
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: MIT License",
			"Operating System :: OS Independent",
		),
		keywords ='genome variants genes transcripts',
		install_requires = requirements,
		zip_safe = False
)

