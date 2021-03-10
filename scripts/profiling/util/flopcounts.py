import sys, os

def usage():
  print('Usage:\n'
      '\tpython flopcounts.py file linestart lineend')

def main() -> int:
  args = sys.argv[1:]
  if len(args) != 3:
    usage()
    return 1
  f, start, end = args
  start, end = int(start), int(end)
  with open(f,'r') as fp:
    lines = fp.readlines()[start-1:end]

  ops = 0
  s = ''.join(lines)
  slist = ["cost_", "cost","->","/**", "**/","/*","*/","**", "-(-", "-(", "(-"] # Order
  for i in slist:
    s = s.replace(i,'')
  s = s.replace('++','')
  s = s.replace('--','')
  sin = s.count('sin')
  Sin = s.count('Sin')
  cos = s.count('cos')
  Cos = s.count('Cos')
  tan = s.count('tan')
  Tan = s.count('Tan')
  sqrt = s.count('sqrt')
  Sqrt = s.count('Sqrt')
  for op in '+-*/':
    ops += s.count(op)

  # import pdb;pdb.set_trace()
  print(f'sin: {sin}\ncos: {cos}\ntan: {tan}\nsqrt: {sqrt}\nSin: {Sin}\nCos: {Cos}\nTan: {Tan}\nSqrt: {Sqrt}\nflps: {ops}')
  return 0

if __name__ == '__main__':
  sys.exit(main())
