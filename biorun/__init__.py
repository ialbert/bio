VERSION = "0.0.1"

import jinja2


def template(text, ctx={}):
    """
    Renders a template with a context.
    """
    tmpl = jinja2.Template(text, trim_blocks=True, lstrip_blocks=True, autoescape=False, undefined=jinja2.StrictUndefined)
    text = tmpl.render(ctx)
    return text
