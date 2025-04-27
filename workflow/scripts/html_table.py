import pandas as pd

url = snakemake.params['url']
attrs = snakemake.params['attrs']
csv_file = str(snakemake.output)

tables = pd.read_html(url, attrs = attrs)
pd.concat(tables).to_csv(csv_file)
