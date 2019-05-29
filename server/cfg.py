import json
import os


BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class CfgError(Exception):
    pass


def get_param(setting):
    try:
        return cfg[setting]

    except KeyError:
        error_msg = "variable {0} wasn't found in config".format(setting)
        raise CfgError(error_msg)


with open(os.path.join(BASE_DIR, "cfg.json")) as f:
    cfg = json.loads(f.read())
