from box.misc import AndersonMixer
from box.pulay import PulayMixer

mixers = {'pulay':PulayMixer, 'anderson':AndersonMixer}

def BuildMixer(params):
    if params == None:
        return mixers['anderson']()
    elif type(params) == str:
        name = params
        return mixers[name.lower()]()
    elif type(params) == dict and 'name' in params:
        p = params.copy()
        name = p['name'].lower()
        p.pop('name')
        return mixers[name](**p)
    else:
        raise Exception('You must provide the name of the mixer.')

