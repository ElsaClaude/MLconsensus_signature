from snakemake.utils import validate
import pandas as pd

##### load config file #####

configfile: "../config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")