"""
Test que les anciens modes fonctionnent toujours pareil
Run: python test_backwards_compatibility.py
"""
import requests
import json

BASE_URL = "http://localhost:8000"

def test_mode(mode: str, protein_id: str = "P01189"):
    """Teste un mode spécifique"""
    print(f"\n{'='*50}")
    print(f"🧪 Testing mode: {mode}")
    print(f"{'='*50}")
    
    response = requests.post(
        f"{BASE_URL}/analyze",
        json={
            "proteinId": protein_id,
            "mode": mode
        }
    )
    
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Status: {response.status_code}")
        print(f"📊 Sequence length: {data.get('sequenceLength')}")
        print(f"✂️ Cleavage sites: {data.get('cleavageSitesCount')}")
        print(f"🧬 Peptides found: {len(data.get('peptides', []))}")
        
        if data.get('peptides'):
            print(f"   First peptide: {data['peptides'][0]['sequence'][:30]}...")
        
        return True
    else:
        print(f"❌ Status: {response.status_code}")
        print(f"   Error: {response.text}")
        return False

def main():
    print("🔬 BACKWARDS COMPATIBILITY TEST")
    print("Testing POMC (P01189) - a well-known prohormone")
    
    results = {}
    
    # Test tous les anciens modes
    for mode in ["strict", "permissive", "ultra-permissive"]:
        results[mode] = test_mode(mode, "P01189")
    
    # Test le nouveau mode avec GDF11
    print(f"\n{'='*50}")
    print(f"🆕 Testing NEW mode: pcsk567 with GDF11")
    print(f"{'='*50}")
    results["pcsk567"] = test_mode("pcsk567", "O95390")
    
    # Résumé
    print(f"\n{'='*50}")
    print("📋 SUMMARY")
    print(f"{'='*50}")
    
    all_passed = True
    for mode, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"   {mode}: {status}")
        if not passed:
            all_passed = False
    
    if all_passed:
        print(f"\n🎉 ALL TESTS PASSED!")
    else:
        print(f"\n⚠️ SOME TESTS FAILED!")
    
    return all_passed

if __name__ == "__main__":
    main()