import re, json

prottest_file = snakemake.input['prottest']
output_file = str(snakemake.output)

criterion = snakemake.params['criterion']

def parse_model(line, file):
    model_name = line[40:].strip()
    model = { 'model': model_name.split('+')[0] }
    for line in file:
        if len(line) < 40: break
        key = line[:37].strip(' .')
        value = line[40:].strip()
        if key.startswith('gamma shape'):
            model['gamma'] = value
        elif key == 'proportion of invariable sites':
            model['pinvar'] = value
        elif key == 'aminoacid frequencies':
            model['freqs'] = 'observed' if value.startswith('observed') else value
    return model_name, model

def parse_prottest(model_file, criterion):
    models = {}
    with open(model_file) as file:
        try:
            while True:
                line = next(file)
                if line.startswith('Model.'):
                    model_name, model = parse_model(line, file)
                    models[model_name] = model
                elif line.startswith('Best model according to'):
                    sel_criterion, model_name = line[24:].split(':')
                    model_name = model_name.strip()
                    return models[model_name]
        except EOFError:
            return None

model = parse_prottest(prottest_file, criterion)
assert model, "Something went wrong with the parsing of the Prottest3 output"

with open(output_file, 'w') as file:    
    json.dump(model, file)
