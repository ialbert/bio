"""
Command-line interface for the bio package.
"""

import click
import sys
from biorun import search_lib

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
    """
    bio: making bioinformatics fun again
    """
    # If no subcommand was provided, show help
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        sys.exit(0)

def auto_help(ctx):
    """
    Show command help if no arguments were provided
    """
    if (len(sys.argv) == 2):
        click.echo(ctx.get_help())
        sys.exit(0)

@cli.command()
@click.pass_context
def fetch(ctx):
    """
    Fetch biological data.
    
    This command allows you to fetch biological data from various sources.
    """
    # Print help if no arguments were provided
    auto_help(ctx)

    click.echo("Fetch command - implementation coming soon")


@cli.command()
@click.pass_context
@click.argument('words', nargs=-1)
def search(ctx, words):
    """
    Search biological data.
    
    This command allows you to search through biological data.
    """
    
    # Print help if no arguments were provided
    auto_help(ctx)
    
    for word in words:
        type = search_lib.get_type(word)
        click.echo(f"{word} - {type}")

    click.echo("Search command - implementation coming soon")

if __name__ == '__main__':
    sys.argv.append("search")
    sys.argv.append("AF086833")
    sys.argv.append("PRJNA257197")
    sys.argv.append("SRR14575325")
    sys.argv.append("GCA_000002415")

    cli()