# Install
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n transcriptionary --file requirements.txt python=3.10
conda activate transcriptionary
python setup.py install
```

or

```
pip install -r requirements.txt
python setup.py install
```

# Testing
```
transcriptionary test/test.yaml
```

will create `transcriptionary-example.html`


# Development

`transcriptionary` is a command-line application. To develop, run:
```
python setup.py develop
```
then changes will be present when running `transcriptionary`