"""
Command-line interface for the bio package.
"""

import click
import sys
from biorun import search_lib
from biorun import utils_lib
from biorun import logger

# Initialize the logger.
log = logger.default

class CustomFormat(click.Group):
    def format_usage(self, ctx, formatter):
        """
        Override usage string.
        """
        formatter.write_usage(ctx.command_path, 'COMMAND [ARGS] ...')
    
    def format_help(self, ctx, formatter):
        """
        Override help output.
        """
        # Call the parent's format_help to get the standard help output
        super().format_help(ctx, formatter)
        
        # Add examples section
        with formatter.section('Examples'):
            formatter.write_text('  bio search SRR14575325')
            formatter.write_text('  bio search PRJNA374918')

@click.group(cls=CustomFormat, invoke_without_command=True, context_settings={'help_option_names': ['-h', '--help']})
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
@click.argument('query')
@click.option('--csv', is_flag=True, help='Output in CSV format')
@click.option('--lines', is_flag=True, help='Line-delimited JSON format')
@click.option('-v', '--verbose', is_flag=True, help='Verbose output')
def search(ctx, query, csv, lines, verbose):
    """
    Search biological data.
    
    This command allows you to search through biological data.
    """
        # Print help if no arguments were provided
    auto_help(ctx)


    # Default logging level.
    if verbose:
        log.setLevel(logger.DEBUG)

    results = search_lib.search(query)
    if csv:
        utils_lib.as_csv(results)
    else:
        results = map(search_lib.update_ftp_fields, results)
        utils_lib.as_json(results, lines=lines)
     


if __name__ == '__main__':
    sys.argv.append("search")
    #sys.argv.append("--csv")
    #sys.argv.append("--lines")
    #sys.argv.append("-v")
    #sys.argv.append("AF086833")
    #sys.argv.append("PRJNA257197")
    sys.argv.append("SRR14575325")
    #sys.argv.append("GCA_000002415")

    cli()