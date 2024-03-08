import os
import pty
import subprocess
import sys

def main():
    """Executes the 'swiftest_driver' binary located in the package root, passing along any command-line arguments, and streams the output to the terminal in real-time, handling progress bars correctly by using a pseudo-terminal."""

    # Determine the path to the binary relative to this script
    package_root = os.path.dirname(os.path.abspath(__file__))
    binary_path = os.path.join(package_root, 'swiftest')

    # sys.argv[1:] contains all the arguments passed to the script, excluding the script name itself
    args = [binary_path] + sys.argv[1:]

    # Use pty to spawn the process and create a pseudo-terminal
    main_fd, subordinate_fd = pty.openpty()

    # Spawn the subprocess
    proc = subprocess.Popen(args, stdin=subordinate_fd, stdout=subordinate_fd, stderr=subordinate_fd, close_fds=True)

    # Close the subordinate file descriptor in the parent process
    os.close(subordinate_fd)

    # Read the output from the main file descriptor
    while True:
        try:
            output = os.read(main_fd, 1024).decode()
            if output:
                sys.stdout.write(output)
                sys.stdout.flush()
            else:
                break
        except OSError:
            break

    # Wait for the subprocess to finish
    proc.wait()

    return

if __name__ == "__main__":
    main()
