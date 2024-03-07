import subprocess
import os
import sys

def run_swiftest_driver():
    """Executes the 'swiftest_driver' binary located in the package root, passing along any command-line arguments."""
    
    # Determine the path to the binary relative to this script
    package_root = os.path.dirname(os.path.abspath(__file__))
    binary_path = os.path.join(package_root, 'swiftest_driver')
    
    # sys.argv[1:] contains all the arguments passed to the script, excluding the script name itself
    args = sys.argv[1:]

    # Execute the binary with additional arguments using subprocess
    result = subprocess.run([binary_path] + args, capture_output=True, text=True)

    # Optionally, print or process the output and errors
    print(result.stdout)
    if result.stderr:
        print(f"Error: {result.stderr}")

    return

if __name__ == "__main__":
    run_swiftest_driver()
