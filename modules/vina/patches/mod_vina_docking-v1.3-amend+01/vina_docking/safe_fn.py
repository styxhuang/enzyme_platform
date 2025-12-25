from hashlib import md5
from time import time

from typing_extensions import Literal

invalid_chars = r'\/:*?"<>|'


def safe_fn(provided_name: str, cata: Literal['Receptor', 'Ligand']) -> str:
    for invalid_char in invalid_chars:
        if invalid_char in provided_name:
            break

    else:
        return provided_name

    hashed = md5(f'{provided_name} {time()}'.encode('utf8'))
    return f'{cata}{hashed.hexdigest().upper()[:6]}'
