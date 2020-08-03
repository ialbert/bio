VERSION = "0.0.1"

from jinja2 import Template


def render(text, ctx={}):
    """
    Renders a template with a context.
    """
    tmpl = Template(text, trim_blocks=True, lstrip_blocks=True, autoescape=False)
    text = tmpl.render(ctx)
    return text
