import os

LOG_LEVEL = os.getenv('BIO') or 'DEBUG'


LOGGING = {

    'version': 1,

    'disable_existing_loggers': True,

    'formatters': {

        'verbose': {
            'format': '%(levelname)s\t%(asctime)s\t%(module)s.%(funcName)s\t%(lineno)s\t%(message)s\t'
        },

        'simple': {
            'format': '%(levelname)s\t%(module)s\t%(message)s'
        },

    },

    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },

    },

    'loggers': {

        'bio': {
            'handlers': ['console'],
            'level': LOG_LEVEL,
        },
        'root': {
             'handlers': ['console'],
             'level': 'DEBUG'
        }

    },
}
