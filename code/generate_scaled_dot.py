#!/usr/bin/env python3
"""
Generate a new dot file from counterexample1_petersen_moebius.dot where:
  - pos= coordinates are kept unchanged (vertex positions stay as-is)
  - edge labels (taillabel/headlabel) are updated to reflect radius-2 scaling

From gp2-1.tex:
  phi = (1 + sqrt(5)) / 2   (golden ratio)
  x   = 2 / 5^(1/4)
  y   = x * phi

The original dot file uses coordinates on the sphere of radius 2*phi.
At radius 2 (scale factor 1/phi), each coordinate c becomes c/phi.

Symbolic simplification rules (dividing by phi), using x/phi = y-x:
  2*phi  ->  2
  phi+1  ->  phi          (since (phi+1)/phi = phi^2/phi = phi)
  phi    ->  1
  1      ->  phi-1        (since 1/phi = phi-1)
  y      ->  x            (since y/phi = x*phi/phi = x)
  x      ->  y-x          (since x/phi = x*(phi-1) = x*phi-x = y-x)

We represent each coordinate as  a + b*phi + c*x + d*y  with integer
coefficients, divide by phi, then format the result canonically.
"""

import re

PHI = "\u03c6"   # φ  (Unicode phi for output)


def parse_expr(s: str):
    """
    Parse a symbolic expression string into integer coefficients (a, b, c, d)
    representing  a + b*phi + c*x + d*y.

    Supported atoms: integers, phi (or φ), x, y.
    Supported structure: sum/difference of terms, each term optionally
    prefixed by an integer coefficient and '*'.
    """
    s = s.strip().replace(PHI, 'phi').replace(' ', '')

    a, b, c, d = 0, 0, 0, 0

    # Ensure a leading sign so every token has an explicit sign
    if s and s[0] not in '+-':
        s = '+' + s

    # Match: sign, optional integer coefficient, optional '*', atom
    pattern = re.compile(
        r'([+\-])'           # sign
        r'(\d+)?'            # optional integer coefficient
        r'(?:\*)?'           # optional '*'
        r'(phi|x|y|\d+)?'   # atom: phi, x, y, or integer
    )

    for m in pattern.finditer(s):
        sign_ch, coeff_str, atom = m.group(1), m.group(2), m.group(3)
        if coeff_str is None and atom is None:
            continue
        sign = 1 if sign_ch == '+' else -1
        coeff = int(coeff_str) if coeff_str else 1

        if atom is None:
            # bare integer coefficient with no atom → constant
            a += sign * coeff
        elif atom == 'phi':
            b += sign * coeff
        elif atom == 'x':
            c += sign * coeff
        elif atom == 'y':
            d += sign * coeff
        else:
            # atom is a digit string → multiply into constant
            a += sign * coeff * int(atom)

    return a, b, c, d


def divide_by_phi(a, b, c, d):
    """
    Divide  a + b*phi + c*x + d*y  by phi, using:
      1/phi   = phi - 1          =>  a  ->  -a + a*phi
      phi/phi = 1                =>  b  ->  b  (constant)
      x/phi   = y - x            =>  c  ->  -c*x + c*y
      y/phi   = x                =>  d  ->  d*x

    Returns new (a', b', c', d').
    """
    new_a = -a + b       # constant part
    new_b = a            # phi part
    new_c = -c + d       # x part
    new_d = c            # y part
    return new_a, new_b, new_c, new_d


def format_expr(a, b, c, d):
    """
    Format  a + b*phi + c*x + d*y  as a canonical string.
    Term order: y, x, phi, constant (most to least "complex").
    Within that order, positive terms are emitted before negative ones
    so the leading term is positive whenever possible.
    Uses φ for phi.
    """
    # Build list of (coeff, symbol) in canonical order
    ordered = [('y', d), ('x', c), (PHI, b), ('', a)]

    pos_terms = []
    neg_terms = []

    for sym, val in ordered:
        if val == 0:
            continue
        if sym == '':
            tok = str(val)
        elif val == 1:
            tok = sym
        elif val == -1:
            tok = '-' + sym
        else:
            tok = f'{val}*{sym}'

        if val > 0:
            pos_terms.append(tok)
        else:
            neg_terms.append(tok)

    terms = pos_terms + neg_terms

    if not terms:
        return '0'

    # Join: first term as-is; subsequent terms with explicit sign
    result = terms[0]
    for t in terms[1:]:
        if t.startswith('-'):
            result += t        # e.g. "φ-1"
        else:
            result += '+' + t
    return result


def simplify_expr(expr: str) -> str:
    """Parse expr, divide by phi, return formatted result."""
    a, b, c, d = parse_expr(expr)
    na, nb, nc, nd = divide_by_phi(a, b, c, d)
    return format_expr(na, nb, nc, nd)


def negate_formatted(expr: str) -> str:
    """Negate a formatted expression string."""
    a, b, c, d = parse_expr(expr)
    return format_expr(-a, -b, -c, -d)


def transform_label(label_str: str) -> str:
    """
    Transform a full label string like "y, -x, 2 / q 5" by dividing each
    coordinate by phi symbolically, keeping the '/ q N' suffix unchanged.
    """
    if ' / q' in label_str:
        coords_part, q_part = label_str.split(' / q', 1)
        q_suffix = ' / q' + q_part
    else:
        coords_part = label_str
        q_suffix = ''

    raw_coords = [c.strip() for c in coords_part.split(',')]

    new_coords = []
    for coord in raw_coords:
        # Handle outer negation: -(expr)
        m = re.fullmatch(r'-\((.+)\)', coord.strip())
        if m:
            inner_result = simplify_expr(m.group(1))
            new_coords.append(negate_formatted(inner_result))
        else:
            new_coords.append(simplify_expr(coord))

    return ', '.join(new_coords) + q_suffix


def process_line(line: str) -> str:
    """Process a single line of the dot file.
    - pos= attributes: left unchanged
    - taillabel= / headlabel= attributes: coordinates divided by phi symbolically
    """
    def replace_label_attr(m):
        attr_name = m.group(1)
        content   = m.group(2)
        new_content = transform_label(content)
        return f'{attr_name}="{new_content}"'

    line = re.sub(r'(taillabel|headlabel)="([^"]+)"', replace_label_attr, line)
    return line


def main():
    input_path  = "dot/counterexample1_petersen_moebius.dot"
    output_path = "dot/counterexample1_petersen_moebius_r2.dot"

    with open(input_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    output_lines = [process_line(line) for line in lines]

    with open(output_path, 'w', encoding='utf-8') as f:
        f.writelines(output_lines)

    print(f"Written to {output_path}")

    # Print a verification table
    print("\nLabel transformations (original -> radius-2):")
    for line in lines:
        m = re.search(r'(taillabel|headlabel)="([^"]+)"', line)
        if m:
            orig = m.group(2)
            new  = transform_label(orig)
            print(f"  {orig!r:55s}  ->  {new!r}")


if __name__ == "__main__":
    main()
