```
conda create -n transcriptionary -c bioconda --file requirements.txt 
conda activate transcriptionary
```

or

```
pip install -r requirements.txt
```

# testing
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

