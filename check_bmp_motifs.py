"""Check why BMP7 and BMP5 are not detected"""
import re
import requests

proteins = {
    "P18075": "BMP7",
    "P22003": "BMP5"
}

pattern_strict = r"R[A-Z](?:K|R)R"  # Notre pattern actuel
pattern_relaxed = r"R[A-Z]{2}R"     # Pattern plus large (RXXR)

for uniprot_id, name in proteins.items():
    print(f"\n{'='*50}")
    print(f"🔍 {name} ({uniprot_id})")
    print(f"{'='*50}")
    
    # Fetch sequence
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    lines = response.text.split('\n')
    sequence = ''.join(lines[1:])
    
    print(f"Sequence length: {len(sequence)} aa")
    
    # Search with strict pattern
    strict_matches = list(re.finditer(pattern_strict, sequence))
    print(f"\nStrict pattern R[A-Z](K|R)R: {len(strict_matches)} matches")
    for m in strict_matches:
        print(f"   - {m.group()} at position {m.start() + 1}")
    
    # Search with relaxed pattern
    relaxed_matches = list(re.finditer(pattern_relaxed, sequence))
    print(f"\nRelaxed pattern R[A-Z][A-Z]R: {len(relaxed_matches)} matches")
    for m in relaxed_matches:
        print(f"   - {m.group()} at position {m.start() + 1}")
    
    # Show last 150 aa (where mature form usually starts)
    print(f"\nLast 150 aa (C-terminal region):")
    print(f"   ...{sequence[-150:]}")