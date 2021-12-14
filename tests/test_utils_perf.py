import datamol as dm


def test_watch_duration():
    def fn(n):
        for i in range(n):
            print(i)

    with dm.utils.perf.watch_duration(log=True) as w:
        fn(5)

    assert isinstance(w.duration, float)
