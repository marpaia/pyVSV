"""Command-line entry point for VSV."""

from __future__ import annotations

from collections.abc import Sequence

import click

from vsv import __version__
from vsv.cli.design import antigenome, check, designs, precedents, show
from vsv.cli.genome import annotate, extract, genes
from vsv.errors import InputError


@click.group(invoke_without_command=True)
@click.version_option(__version__, message="%(prog)s %(version)s")
@click.pass_context
def cli(ctx: click.Context) -> None:
    """Computational genetics tools for vesicular stomatitis virus research."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


cli.add_command(annotate)
cli.add_command(genes)
cli.add_command(extract)
cli.add_command(designs)
cli.add_command(check)
cli.add_command(show)
cli.add_command(precedents)
cli.add_command(antigenome)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI and return its exit status."""
    try:
        result = cli.main(
            args=None if argv is None else list(argv),
            prog_name="vsv",
            standalone_mode=False,
        )
    except InputError as error:
        click.echo(f"vsv: error: {error}", err=True)
        return 1
    except click.ClickException as error:
        error.show()
        return error.exit_code
    return result if isinstance(result, int) else 0
