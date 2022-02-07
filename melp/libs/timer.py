import time


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError("Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self, verbosity=False):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        if verbosity is False:
            print(f"Elapsed time: {elapsed_time:0.4f} seconds")

    def get_elapsed_time(self, verbosity=False):
        """return the elapsed time"""
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        if verbosity is False:
            print(f"Elapsed time: {elapsed_time:0.4f} seconds")
        return elapsed_time

    def print(self):
        self.get_elapsed_time(verbosity=False)
