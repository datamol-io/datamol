import multiprocessing


DATAMOL_MAX_IMPORT_DURATION = 2  # in seconds


def _import_datamol_fn(queue):
    import datetime

    start = datetime.datetime.now()
    import datamol

    duration = datetime.datetime.now() - start
    duration = duration.total_seconds()
    queue.put(duration)


def test_datamol_import():

    context = multiprocessing.get_context(method="spawn")

    queue = context.Queue()
    p = context.Process(target=_import_datamol_fn, args=(queue,))
    p.start()
    duration = queue.get()
    p.join()
    p.close()

    print(duration)

    # Check datamol wasn't too fast to load
    # That should never happen and it's a simple sanity check
    # whether the `spawn` process worked correctly.
    assert duration > 1e-2, "datamol loaded too fast, something is wrong with the test."

    # Check datamol wasn't too long to load
    assert duration < DATAMOL_MAX_IMPORT_DURATION, "datamol tooks too long to load."
