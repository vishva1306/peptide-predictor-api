"""
Test PCSK5/6/7 detection on GDF11
Run: python test_pcsk567.py
"""
import re

# ==================== GDF11 SEQUENCE (O95390) ====================
# Source: UniProt O95390
GDF11_SEQUENCE = """
MVLAAPLLLGFLLLALEL RPRGEAAEGP AAAAAAAAAA AGVGGERSSR PAPSVAPEPD
GCPVCVWRQH SRELRLESI KSQILSKLRL KEAPNISREV VKQLLPKAPP LQQILDLHDF
QGDALQPEDF LEEDEYHAT TETVISMAQE TDPAVQTDGS PLCCHFHFSP KVMFTKVLKA
QLWVYLRPVP RPATVYLQIL RLKPLTGEGT AGGGGGGRRH IRIRSLKIEL HSRSGHWQSI
DFKQVLHSWF RQPQSNWGIE INAFDPSGTD LAVTSLGPGA EGLHPFMELR VLENTKRSRR
NLGLDCDEH SSESRCCRYP LTVDFEAFGW DWIIAPKRYK ANYCSGQCEY MFMQKYPTHT
HLVQQANPRG SAGPCCTPTK MSPINMLYFN DKQQIIYGKI PGMVVDRCGC S
""".replace("\n", "").replace(" ", "")

print("=" * 60)
print("🔬 TEST PCSK5/6/7 DETECTION ON GDF11 (O95390)")
print("=" * 60)

print(f"\n📊 Sequence length: {len(GDF11_SEQUENCE)} aa")
print(f"📋 First 50 aa: {GDF11_SEQUENCE[:50]}...")
print(f"📋 Last 50 aa: ...{GDF11_SEQUENCE[-50:]}")

# ==================== PATTERN PCSK5/6/7 ====================
# R-X-(K/R)-R where X = any amino acid
pattern = r"R[A-Z](?:K|R)R"

print(f"\n🔍 Pattern: {pattern}")
print("   R = Arginine")
print("   [A-Z] = Any amino acid")
print("   (?:K|R) = Lysine or Arginine")
print("   R = Arginine")

# ==================== FIND MATCHES ====================
matches = list(re.finditer(pattern, GDF11_SEQUENCE))

print(f"\n✅ Found {len(matches)} PCSK5/6/7 site(s):")

for i, match in enumerate(matches, 1):
    motif = match.group()
    position = match.start()
    cleavage_after = match.end()
    
    print(f"\n   Site {i}:")
    print(f"   - Motif: {motif}")
    print(f"   - Position: {position + 1} (1-indexed)")
    print(f"   - Cleavage after position: {cleavage_after}")
    
    # Context around the site
    context_start = max(0, position - 10)
    context_end = min(len(GDF11_SEQUENCE), position + 15)
    context = GDF11_SEQUENCE[context_start:context_end]
    marker_pos = position - context_start
    
    print(f"   - Context: ...{context[:marker_pos]}[{context[marker_pos:marker_pos+4]}]{context[marker_pos+4:]}...")
    
    # Extract mature form (after cleavage)
    mature_seq = GDF11_SEQUENCE[cleavage_after:]
    print(f"   - Mature form: {len(mature_seq)} aa")
    print(f"   - First 50 aa of mature: {mature_seq[:50]}...")

# ==================== VALIDATION ====================
print("\n" + "=" * 60)
print("🧪 VALIDATION")
print("=" * 60)

# Expected: RSRR should be found (from literature)
expected_motifs = ["RSRR", "RKRR", "RSKR"]
found_motifs = [m.group() for m in matches]

print(f"\n📋 Expected motifs (from literature): {expected_motifs}")
print(f"📋 Found motifs: {found_motifs}")

# Check for RSRR specifically
if "RSRR" in found_motifs:
    print("\n✅ SUCCESS: Found RSRR motif (expected for GDF11)")
    
    # Find the RSRR match
    for match in matches:
        if match.group() == "RSRR":
            mature_start = match.end()
            mature_seq = GDF11_SEQUENCE[mature_start:]
            
            print(f"\n📊 MATURE GDF11 DETAILS:")
            print(f"   - Start position: {mature_start + 1}")
            print(f"   - Length: {len(mature_seq)} aa")
            print(f"   - Expected length: ~109 aa (from literature)")
            
            # Check if length is approximately correct
            if 100 <= len(mature_seq) <= 120:
                print(f"   ✅ Length is correct!")
            else:
                print(f"   ⚠️ Length differs from expected")
            
            print(f"\n   MATURE SEQUENCE:")
            print(f"   {mature_seq}")
else:
    print("\n❌ FAILED: RSRR motif not found")
    print("   Check if the sequence is correct")

# ==================== ADDITIONAL CHECKS ====================
print("\n" + "=" * 60)
print("🔬 ADDITIONAL CHECKS")
print("=" * 60)

# Check for signal peptide
signal_peptide = GDF11_SEQUENCE[:18]
print(f"\n📍 Signal peptide (first 18 aa): {signal_peptide}")

# Check for prodomain
for match in matches:
    if match.group() == "RSRR":
        prodomain = GDF11_SEQUENCE[18:match.start()]
        print(f"📍 Prodomain (after signal, before RSRR): {len(prodomain)} aa")
        print(f"   {prodomain[:50]}...")
        break

print("\n" + "=" * 60)
print("✅ TEST COMPLETED")
print("=" * 60)