"""Allow ``python -m fpcaller`` as an alternative to the ``fpcaller`` command."""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
