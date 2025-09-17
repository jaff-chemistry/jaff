#!/usr/bin/env python3
# ABOUTME: Script to extract and test Python code blocks from markdown documentation files
# ABOUTME: Used by GitHub Actions to ensure all documentation examples are runnable

import argparse
import ast
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple


def extract_python_blocks(markdown_content: str, filename: str) -> List[Tuple[str, int]]:
    """Extract Python code blocks from markdown content.
    
    Returns list of (code, line_number) tuples.
    """
    # Pattern to match Python code blocks
    pattern = r'```python\n(.*?)\n```'
    blocks = []
    
    for match in re.finditer(pattern, markdown_content, re.DOTALL):
        code = match.group(1)
        # Calculate line number where code block starts
        line_num = markdown_content[:match.start()].count('\n') + 2
        blocks.append((code, line_num))
    
    return blocks


def filter_testable_code(code: str) -> Tuple[bool, str]:
    """Determine if code block should be tested and clean it.
    
    Returns (should_test, cleaned_code).
    """
    lines = code.strip().split('\n')
    
    # Skip pure shell/bash commands (even if in python block)
    if all(line.strip().startswith('#') or line.strip() == '' or 
           line.strip().startswith('jaff ') for line in lines):
        return False, code
    
    # Check if it's a code snippet without complete context
    has_ellipsis = any('...' in line for line in lines)
    if has_ellipsis:
        return False, code
        
    # Skip code blocks with control flow statements that require surrounding context
    control_flow_keywords = ['continue', 'break']
    for line in lines:
        stripped = line.strip()
        if any(stripped == keyword or stripped.startswith(keyword + ' ') 
               for keyword in control_flow_keywords):
            return False, code
        
    # Filter out matplotlib.pyplot.show() calls to avoid GUI issues in CI
    # and fix common documentation issues
    cleaned_lines = []
    for line in lines:
        if 'plt.show()' in line:
            cleaned_lines.append(line.replace('plt.show()', 'plt.close()'))
        elif 'rates.shape' in line and 'rates = network.get_table(' in ''.join(lines[:lines.index(line)]):
            # If we're accessing .shape on rates after get_table, add temp variable
            cleaned_lines.append(line)
        else:
            cleaned_lines.append(line)
    
    cleaned_code = '\n'.join(cleaned_lines)
    return True, cleaned_code


def test_file_examples(md_file: Path, test_dir: Path, verbose: bool = False) -> Tuple[int, int, List[str]]:
    """Test all Python examples in a markdown file together.
    
    Returns (total_blocks, failed_count, error_messages).
    """
    content = md_file.read_text()
    blocks = extract_python_blocks(content, str(md_file))
    
    if not blocks:
        return 0, 0, []
    
    # Collect all testable code blocks
    all_code_parts = []
    block_info = []  # Track which lines correspond to which blocks
    
    for i, (code, line_num) in enumerate(blocks):
        should_test, cleaned_code = filter_testable_code(code)
        if should_test:
            all_code_parts.append(f"# Block {i+1} from line {line_num}")
            all_code_parts.append(cleaned_code)
            all_code_parts.append("")  # Empty line between blocks
            block_info.append((i+1, line_num))
    
    if not all_code_parts:
        if verbose:
            print(f"  No testable Python blocks found")
        return len(blocks), 0, []
    
    # Combine all code blocks
    combined_code = '\n'.join(all_code_parts)
    
    # Create test file with setup
    test_file = test_dir / "test_combined.py"
    
    # Get absolute path to project networks directory
    project_root = Path('.').resolve()
    source_networks_dir = project_root / 'networks'
    
    setup_code = f"""
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

# Suppress matplotlib display backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Create networks directory and copy real network files if they exist
import shutil
os.makedirs('networks', exist_ok=True)

# Copy any existing network files to temp directory
source_networks_dir = r'{source_networks_dir}'
if os.path.exists(source_networks_dir):
    for file in os.listdir(source_networks_dir):
        src_path = os.path.join(source_networks_dir, file)
        dest_path = os.path.join('networks', file)
        if os.path.isfile(src_path):
            shutil.copy2(src_path, dest_path)

# Create symbolic links to real network files that exist
network_files_to_link = [
    ('networks/gas_reactions_kida.uva.2024.in', 'gas_reactions_kida.uva.2024.in'),
    ('networks/react_popsicle_semenov', 'react_popsicle_semenov'),
    ('networks/react_COthin', 'react_COthin')
]

for target_path, source_file in network_files_to_link:
    source_path = os.path.join(source_networks_dir, source_file)
    if os.path.exists(source_path) and not os.path.exists(target_path):
        # Create symbolic link only if source exists
        os.symlink(source_path, target_path)

# Create output files that might be generated
os.makedirs('output', exist_ok=True)

"""
    
    full_code = setup_code + combined_code
    
    with open(test_file, 'w') as f:
        f.write(full_code)
    
    # Run the combined code
    try:
        result = subprocess.run(
            ['python', str(test_file)],
            capture_output=True,
            text=True,
            timeout=900,  # Longer timeout for combined execution
            cwd=str(test_dir)
        )
        
        if result.returncode != 0:
            error_msg = f"Failed to run examples from {md_file}\n"
            error_msg += f"Exit code: {result.returncode}\n"
            if result.stdout:
                error_msg += f"STDOUT:\n{result.stdout}\n"
            if result.stderr:
                error_msg += f"STDERR:\n{result.stderr}\n"
            return len(blocks), 1, [error_msg]
        
        if verbose:
            print(f"  ✓ All {len(block_info)} testable blocks passed")
        
        return len(blocks), 0, []
        
    except subprocess.TimeoutExpired:
        return len(blocks), 1, [f"Execution timed out for {md_file} (900s limit)"]
    except Exception as e:
        return len(blocks), 1, [f"Unexpected error testing {md_file}: {str(e)}"]


def test_code_block(code: str, test_dir: Path) -> Tuple[bool, str]:
    """Test a Python code block by running it.
    
    Returns (success, error_message).
    """
    # Create a test file
    test_file = test_dir / "test_example.py"
    
    # Add necessary imports and setup
    setup_code = """
import sys
import os
sys.path.insert(0, os.path.abspath('.'))

# Suppress matplotlib display backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Create networks directory if referenced
if 'networks/' in '''{}''':
    os.makedirs('networks', exist_ok=True)
    # Create dummy network files for testing
    with open('networks/react_COthin', 'w') as f:
        f.write('# Dummy network file\\n')
    with open('networks/react_popsicle_semenov', 'w') as f:
        f.write('# Dummy network file\\n')
    with open('networks/gas_reactions_kida.uva.2024.in', 'w') as f:
        f.write('# Dummy KIDA network\\n')
    with open('networks/prizmo_network.dat', 'w') as f:
        f.write('# Dummy PRIZMO network\\n')
    with open('networks/krome_network.dat', 'w') as f:
        f.write('# Dummy KROME network\\n')

""".format(code)
    
    full_code = setup_code + code
    
    with open(test_file, 'w') as f:
        f.write(full_code)
    
    # Run the code
    try:
        result = subprocess.run(
            ['python', str(test_file)],
            capture_output=True,
            text=True,
            timeout=900,
            cwd=str(test_dir)
        )
        
        if result.returncode != 0:
            return False, f"Exit code {result.returncode}\\nSTDOUT:\\n{result.stdout}\\nSTDERR:\\n{result.stderr}"
        
        return True, ""
        
    except subprocess.TimeoutExpired:
        return False, "Code execution timed out (900s limit)"
    except Exception as e:
        return False, f"Unexpected error: {str(e)}"


def main():
    parser = argparse.ArgumentParser(description='Test Python examples in documentation')
    parser.add_argument('files', nargs='*', help='Markdown files to test (default: all docs/*.md)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show all tested code blocks')
    parser.add_argument('--fail-fast', action='store_true', help='Stop on first failure')
    args = parser.parse_args()
    
    # Find markdown files
    if args.files:
        md_files = args.files
    else:
        docs_dir = Path('docs')
        if not docs_dir.exists():
            print(f"Error: docs directory not found", file=sys.stderr)
            sys.exit(1)
        md_files = list(docs_dir.glob('**/*.md'))
    
    total_blocks = 0
    total_files = 0
    failed_files = 0
    all_errors = []
    
    # Create temporary directory for tests
    with tempfile.TemporaryDirectory() as temp_dir:
        test_dir = Path(temp_dir)
        
        for md_file in md_files:
            md_path = Path(md_file)
            if not md_path.exists():
                print(f"Warning: {md_file} not found", file=sys.stderr)
                continue
            
            total_files += 1
            print(f"\nTesting examples in {md_file}...")
            
            blocks_count, failed_count, errors = test_file_examples(md_path, test_dir, args.verbose)
            total_blocks += blocks_count
            
            if failed_count > 0:
                failed_files += 1
                all_errors.extend(errors)
                
                for error in errors:
                    print(f"  ❌ FAILED")
                    print(f"    {error}")
                
                if args.fail_fast:
                    sys.exit(1)
            elif blocks_count == 0:
                print(f"  No Python code blocks found")
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Total files tested: {total_files}")
    print(f"Total code blocks found: {total_blocks}")
    print(f"Files with failures: {failed_files}")
    
    if failed_files > 0:
        print(f"\nFailed files:")
        for error in all_errors:
            print(f"  - {error.split(chr(10))[0]}")
        sys.exit(1)
    else:
        print("\nAll tests passed! ✓")


if __name__ == '__main__':
    main()