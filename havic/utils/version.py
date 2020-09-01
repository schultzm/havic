#!/usr/bin/env python3

from .. import (__version__, __version_date__, __author__, __author_email__,
                __github_username__, __download_url__, )


class Version:
    """Print the version to stdout."""

    def __init__(self):
        print("Version:", __version__)
        print("Version date:", __version_date__)
        print("Author:", __author__)
        print("Author email:", __author_email__)
        print("Github username:", __github_username__)
        print("Download url:", __download_url__)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
