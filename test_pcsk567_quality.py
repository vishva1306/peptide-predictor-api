"""
Test PCSK5/6/7 quality - Known substrates from literature
Run: python test_pcsk567_quality.py
"""
import requests
import json

BASE_URL = "http://localhost:8000"

# ==================== KNOWN PCSK5/6/7 SUBSTRATES ====================
# Source: Table 2 from the paper mentioned by your colleague
KNOWN_SUBSTRATES = [
    {
        "id": "O95390",
        "name": "GDF11 (Growth/differentiation factor 11)",
        "expected_motif": "RSRR",
        "expected_mature_length": (100, 120),  # ~109 aa
        "expected_mature_start": "NLGLD",  # First 5 aa of mature form
    },
    {
        "id": "P01137",
        "name": "TGFB1 (Transforming growth factor beta-1)",
        "expected_motif": "RHRR",  # R-X-R-R pattern
        "expected_mature_length": (100, 130),
        "expected_mature_start": None,  # Will check if found
    },
    {
        "id": "P18075",
        "name": "BMP7 (Bone morphogenetic protein 7)",
        "expected_motif": None,  # Will check what we find
        "expected_mature_length": (100, 150),
        "expected_mature_start": None,
    },
    {
        "id": "P12643",
        "name": "BMP2 (Bone morphogenetic protein 2)",
        "expected_motif": None,
        "expected_mature_length": (100, 130),
        "expected_mature_start": None,
    },
    {
        "id": "P22003",
        "name": "BMP5 (Bone morphogenetic protein 5)",
        "expected_motif": None,
        "expected_mature_length": (100, 150),
        "expected_mature_start": None,
    },
]

def test_substrate(substrate: dict) -> dict:
    """Test a known PCSK5/6/7 substrate"""
    print(f"\n{'='*60}")
    print(f"🧪 Testing: {substrate['name']}")
    print(f"   UniProt ID: {substrate['id']}")
    print(f"{'='*60}")
    
    result = {
        "id": substrate["id"],
        "name": substrate["name"],
        "status": "UNKNOWN",
        "details": {}
    }
    
    try:
        response = requests.post(
            f"{BASE_URL}/analyze",
            json={
                "proteinId": substrate["id"],
                "mode": "pcsk567"
            },
            timeout=30
        )
        
        if response.status_code == 200:
            data = response.json()
            
            print(f"\n📊 Results:")
            print(f"   Sequence length: {data.get('sequenceLength')} aa")
            print(f"   Cleavage sites found: {data.get('cleavageSitesCount')}")
            print(f"   Peptides extracted: {len(data.get('peptides', []))}")
            
            result["details"]["sequenceLength"] = data.get("sequenceLength")
            result["details"]["cleavageSites"] = data.get("cleavageSitesCount")
            result["details"]["peptides"] = len(data.get("peptides", []))
            
            # Check cleavage sites
            sites = data.get("cleavageSites", [])
            if sites:
                print(f"\n✂️ Cleavage sites:")
                for site in sites:
                    print(f"   - {site['motif']} at position {site['position']}")
                
                # Check expected motif
                if substrate["expected_motif"]:
                    found_motifs = [s["motif"] for s in sites]
                    if substrate["expected_motif"] in found_motifs:
                        print(f"   ✅ Expected motif {substrate['expected_motif']} FOUND!")
                        result["details"]["motif_match"] = True
                    else:
                        print(f"   ⚠️ Expected {substrate['expected_motif']}, found {found_motifs}")
                        result["details"]["motif_match"] = False
            
            # Check peptides
            peptides = data.get("peptides", [])
            mature_forms = [p for p in peptides if p.get("peptideType") == "mature_form"]
            
            if mature_forms:
                print(f"\n🧬 Mature form(s):")
                for mf in mature_forms:
                    print(f"   - Length: {mf['length']} aa")
                    print(f"   - Start: {mf['sequence'][:20]}...")
                    print(f"   - Position: {mf['start']} → {mf['end']}")
                    
                    # Check expected length
                    if substrate["expected_mature_length"]:
                        min_len, max_len = substrate["expected_mature_length"]
                        if min_len <= mf["length"] <= max_len:
                            print(f"   ✅ Length {mf['length']} is within expected range ({min_len}-{max_len})")
                            result["details"]["length_match"] = True
                        else:
                            print(f"   ⚠️ Length {mf['length']} outside expected range ({min_len}-{max_len})")
                            result["details"]["length_match"] = False
                    
                    # Check expected start
                    if substrate["expected_mature_start"]:
                        if mf["sequence"].startswith(substrate["expected_mature_start"]):
                            print(f"   ✅ Starts with expected '{substrate['expected_mature_start']}'")
                            result["details"]["start_match"] = True
                        else:
                            print(f"   ⚠️ Expected start '{substrate['expected_mature_start']}', got '{mf['sequence'][:5]}'")
                            result["details"]["start_match"] = False
                
                result["status"] = "SUCCESS"
            else:
                print(f"\n⚠️ No mature form found!")
                result["status"] = "NO_MATURE_FORM"
        
        elif response.status_code == 404:
            print(f"❌ Protein not found (might not be secreted)")
            result["status"] = "NOT_FOUND"
        else:
            print(f"❌ Error: {response.status_code}")
            print(f"   {response.text[:200]}")
            result["status"] = "ERROR"
    
    except Exception as e:
        print(f"❌ Exception: {e}")
        result["status"] = "EXCEPTION"
    
    return result

def main():
    print("="*60)
    print("🔬 PCSK5/6/7 QUALITY TEST")
    print("   Testing known substrates from literature")
    print("="*60)
    
    results = []
    
    for substrate in KNOWN_SUBSTRATES:
        result = test_substrate(substrate)
        results.append(result)
    
    # Summary
    print("\n" + "="*60)
    print("📋 SUMMARY")
    print("="*60)
    
    success_count = 0
    for r in results:
        status_emoji = {
            "SUCCESS": "✅",
            "NO_MATURE_FORM": "⚠️",
            "NOT_FOUND": "❌",
            "ERROR": "❌",
            "EXCEPTION": "❌",
            "UNKNOWN": "❓"
        }.get(r["status"], "❓")
        
        print(f"   {status_emoji} {r['id']} ({r['name'].split('(')[0].strip()}): {r['status']}")
        
        if r["status"] == "SUCCESS":
            success_count += 1
            details = r.get("details", {})
            if details.get("motif_match") is False:
                print(f"      └─ ⚠️ Motif mismatch")
            if details.get("length_match") is False:
                print(f"      └─ ⚠️ Length outside expected range")
    
    print(f"\n🎯 Success rate: {success_count}/{len(results)} ({100*success_count/len(results):.0f}%)")
    
    if success_count == len(results):
        print("\n🎉 ALL TESTS PASSED!")
    elif success_count > 0:
        print(f"\n⚠️ {len(results) - success_count} test(s) need attention")
    else:
        print("\n❌ All tests failed - check the algorithm")

if __name__ == "__main__":
    main()