from cobra.io.web.load import load_model
from cobra.io import write_sbml_model

if __name__=='__main__':
    model = load_model('iJO1366')
    write_sbml_model(model, 'testdata/iJO1366.xml')