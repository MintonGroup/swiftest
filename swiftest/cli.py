import os
import sys

def main_caf():
    package_root = os.path.dirname(os.path.abspath(__file__))
    binary_path = os.path.join(package_root, "swiftest_caf")
    if os.path.exists(binary_path):
        main("swiftest_caf")
    else:
        print("This version of swiftest has not been compiled with coarray support. The standard version of swiftest will be executed instead.")
        main()
    return

def main(binary_name="swiftest"):
    """
    Executes the 'swiftest' binary located in the package root, passing along any command-line arguments, and streams the output to the terminal in real-time, handling progress bars correctly by using a pseudo-terminal.
    
    Parameters
    ----------
    binary_name : str, optional
        The name of the binary to execute. Default is 'swiftest'.
    """

    # Determine the path to the binary relative to this script
    package_root = os.path.dirname(os.path.abspath(__file__))
    binary_path = os.path.join(package_root, binary_name)

    # Command-line arguments
    args = [binary_path] + sys.argv[1:]

    # Check if we're on Windows
    if os.name == 'nt':
        # Use pywinpty on Windows
        from winpty import PtyProcess
        import threading

        # Create a new pseudo-terminal
        cmd = f"{' '.join(args)}.exe"
        proc = PtyProcess.spawn(cmd)

        def read_output():
            while proc.isalive():
                data = proc.read()
                sys.stdout.write(data)
                sys.stdout.flush()

        t = threading.Thread(target=read_output)
        t.start()

        # Wait for the process to finish
        proc.wait()
        t.join()
    else:
        # Use pty on POSIX systems
        import pty
        import select

        pid, fd = pty.fork()
        if pid == 0:
            # Child process
            os.execv(binary_path, args)
        else:
            # Parent process
            try:
                while True:
                    r, _, _ = select.select([fd], [], [])
                    if fd in r:
                        output = os.read(fd, 1024)
                        if not output:
                            break
                        sys.stdout.buffer.write(output)
                        sys.stdout.flush()
            except OSError:
                pass
            os.waitpid(pid, 0)

if __name__ == "__main__":
    main()