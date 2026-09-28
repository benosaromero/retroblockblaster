import argparse
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", default="dball")
    parser.add_argument("run")
    args = parser.parse_args()

    subprocess.check_call([sys.executable, "retro.py", args.run])
    subprocess.check_call([sys.executable, "bbs.py", "-db", args.db, args.run])
    subprocess.check_call([sys.executable, "make.py", args.run])
    subprocess.check_call([sys.executable, "rank.py", args.run])

if __name__ == "__main__":
    main()
