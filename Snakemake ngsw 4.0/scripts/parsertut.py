import argparse

parser=argparse.ArgumentParser()
parser.add_argument("square", type=int, help="display the square of a given number")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="count")
#parse_args() method actually returns some data from the options specified, in this case, echo.
args=parser.parse_args()
answer= args.square**2
if args.verbosity>=2:
        print("the square of {} equals {}".format(args.square, answer))
elif args.verbosity>=1:
        print("{}^2 == {}".format(args.square, answer))
else:
        print(answer)

##this first example does almost nothing 

##continue the documentation  : https://docs.python.org/2/howto/argparse.html

