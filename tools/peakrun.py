"""Run a command, pass its output through, then report wall time and peak RSS of the child."""
import resource, subprocess, sys, time
t0 = time.time()
rc = subprocess.call(sys.argv[1:])
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
print(f'[peakrun] wall {time.time() - t0:.1f} s  peak RSS {ru.ru_maxrss / 1024:.0f} MB  '
      f'user {ru.ru_utime:.1f} s  sys {ru.ru_stime:.1f} s', flush=True)
sys.exit(rc)
