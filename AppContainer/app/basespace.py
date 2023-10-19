import requests


def basespace_get(uri: str, **params) -> requests.Response:
    return requests.get(uri, headers={'x-access-token': '52d54105977045bb9fc07bcf8ee1d836'}, params=params)