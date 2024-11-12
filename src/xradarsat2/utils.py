import xradarsat2
import logging
import os
import yaml
import re

def get_glob(strlist):
    # from list of str, replace diff by '?'
    def _get_glob(st):
        stglob = ''.join(
            [
                '?' if len(charlist) > 1 else charlist[0]
                for charlist in [list(set(charset)) for charset in zip(*st)]
            ]
        )
        return re.sub(r'\?+', '*', stglob)

    strglob = _get_glob(strlist)
    if strglob.endswith('*'):
        strglob += _get_glob(s[::-1] for s in strlist)[::-1]
        strglob = strglob.replace('**', '*')

    return strglob


def load_config():
    """

    Returns:
        conf: dict
    """
    local_config_path = os.path.join(os.path.dirname(xradarsat2.__file__), 'localxradarsat2-config.yaml')

    if os.path.exists(local_config_path):
        config_path = local_config_path
    else:
        config_path = os.path.join(os.path.dirname(xradarsat2.__file__), 'xradarsat2-config.yaml')

    logging.info('config path: %s', config_path)
    stream = open(config_path, 'r')
    conf = yaml.load(stream, Loader=yaml.CLoader)
    return conf