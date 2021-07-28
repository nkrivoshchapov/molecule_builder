import logging, builtins

def createLogger(classname):
    logger = logging.getLogger(classname)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(name)s:%(levelname)s  %(message)s")

    file_handler = logging.FileHandler("mylog.log")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger

