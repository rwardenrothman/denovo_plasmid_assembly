import requests


def basespace_get(uri: str, **params) -> requests.Response:
    """
        Sends a GET request to the specified URI with additional parameters and an
        access token for authorization. This function is specifically tailored to
        work with APIs like the BaseSpace API.

        Args:
            uri (str): The URI to which the GET request is sent.
            **params: Arbitrary keyword arguments representing additional parameters
                for the GET request.

        Returns:
            requests.Response: The response object containing the server's response
            to the HTTP request.
    """
    return requests.get(uri, headers={'x-access-token': '52d54105977045bb9fc07bcf8ee1d836'}, params=params)