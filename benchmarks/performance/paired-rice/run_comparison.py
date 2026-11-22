import argparse
import ctypes as c
from ctypes import wintypes as w
import json, os, statistics, subprocess, sys, time
from pathlib import Path
class Memory(c.Structure):
    _fields_ = [('cb', w.DWORD), ('faults',w.DWORD)] + [(name,c.c_size_t) for name in
       ('peak_ws','ws','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile')]
get = c.WinDLL('psapi').GetProcessMemoryInfo
get.argtypes = [w.HANDLE,c.POINTER(Memory),w.DWORD]
get.restype = w.BOOL
parser = argparse.ArgumentParser(description='Run the paired BioKt and Biopython rice benchmark on Windows.')
parser.add_argument('dataset', type=Path, help='Directory containing the three rice benchmark inputs')
parser.add_argument('libs', type=Path, help='BioKt installed distribution library directory')
parser.add_argument('output', type=Path, help='Directory for classes and measurement JSON')
parser.add_argument('--deps', type=Path, help='Optional isolated Python dependency directory')
parser.add_argument('--cold-runs', type=int, default=7)
parser.add_argument('--warm-runs', type=int, default=5)
parser.add_argument('--batches', type=int, default=4000)
args = parser.parse_args()
if min(args.cold_runs, args.warm_runs, args.batches) < 1:
    parser.error('Run and batch counts must be positive.')
work = args.output.resolve()
work.mkdir(parents=True, exist_ok=True)
scripts = Path(__file__).resolve().parent
jdk = Path(os.environ['JAVA_HOME'])/'bin'
classpath = str(args.libs.resolve()/'*')
subprocess.run([str(jdk/'javac.exe'), '-cp', classpath, '-d', str(work),
                str(scripts/'TimingBioKt.java')], check=True)
source = str(args.dataset.resolve())
commands = {
    'BioKt': [str(jdk/'java.exe'),'-XX:ActiveProcessorCount=1','-cp',f'{work};{classpath}','TimingBioKt',source],
    'Biopython': [sys.executable,str(scripts/'timing_biopython.py'),str(args.deps.resolve()) if args.deps else '',source],
}
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
def run(cmd,batches):
    memory = Memory()
    memory.cb = c.sizeof(memory)
    peak = 0
    start = time.perf_counter()
    with subprocess.Popen(cmd + [str(batches)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env) as proc:
        while True:
            if get(int(proc._handle),c.byref(memory),memory.cb):
                peak = max(peak,memory.peak_ws)
            if proc.poll() is not None:
                break
            time.sleep(.003)
        elapsed = (time.perf_counter()-start)*1000
        output, error = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(error or output)
    lines = output.splitlines()
    assert lines[0] == '[437255 transcript:Os09t0103500-01 CTC>CGC L16R]', lines
    assert peak
    return {'process_ms':elapsed, 'peak_working_set_mib':peak/2**20,
            'warm_ms_per_50':float(lines[1]) if batches else None, 'call':lines[0]}
result = {name: {'cold':[], 'warm':[]} for name in commands}
for mode, repeats, batches in [('cold',args.cold_runs,0),('warm',args.warm_runs,args.batches)]:
    for i in range(repeats):
        order = list(commands)
        if i%2:
            order.reverse()
        for name in order:
            row = run(commands[name],batches)
            result[name][mode].append(row)
            print(mode, name, i+1, row, flush=True)
(work/'paired-rice-timing.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
for name, metrics in result.items():
    cold, warm = metrics['cold'], metrics['warm']
    for rows, keys in [(cold,['process_ms','peak_working_set_mib']), (warm,['warm_ms_per_50','peak_working_set_mib'])]:
        for key in keys:
            values=[r[key] for r in rows]
            print(name,key,'median',round(statistics.median(values),4),'range',[round(min(values),4),round(max(values),4)])