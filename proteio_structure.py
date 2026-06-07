import time
import ssl
import sys
import socket  # Import socket to override network behavior
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList

# --- FORCE PYTHON TO USE IPv4 ONLY ---
# This bypasses the IPv6 timeout issue shown in your ping test
old_getaddrinfo = socket.getaddrinfo
def new_getaddrinfo(*args, **kwargs):
    responses = old_getaddrinfo(*args, **kwargs)
    # Filter out IPv6 addresses (AF_INET6), keep only IPv4 (AF_INET)
    return [r for r in responses if r[0] == socket.AF_INET]
socket.getaddrinfo = new_getaddrinfo
# --------------------------------------

# Bypass local SSL verification issues
ssl._create_default_https_context = ssl._create_unverified_context

def search_and_download_pdb(sequence, evalue_threshold=1e-5):
    print("🚀 Submitting sequence to online NCBI BLAST (Forced IPv4)...", flush=True)
    print("⚠️  Note: This takes 1 to 5 minutes. Please do not close this window.", flush=True)
    print("Connecting", end="", flush=True)
    
    try:
        result_handle = NCBIWWW.qblast("blastp", "pdb", sequence)
    except Exception as e:
        print(f"\n❌ Network error occurred: {e}")
        return
        
    print("\n✅ Data received! Parsing results...", flush=True)
    
    # Parse the XML results
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
    
    if not blast_record.alignments:
        print("❌ No matching PDB structures found for this sequence.", flush=True)
        return

    print(f"\n✨ Found {len(blast_record.alignments)} matching structures. Top results:\n", flush=True)
    print(f"{'PDB ID':<10}{'Description':<50}{'E-value':<12}{'Similarity (%)':<15}", flush=True)
    print("-" * 90, flush=True)
    
    best_pdb_id = None
    
    for alignment in blast_record.alignments:
        title_parts = alignment.title.split('|')
        pdb_id = title_parts[1].upper() if len(title_parts) > 1 else alignment.hit_id[:4].upper()
        
        for hsp in alignment.hsps:
            if hsp.expect < evalue_threshold:
                identity_percentage = (hsp.identities / hsp.align_length) * 100
                desc = alignment.title.split(' ', 1)[1][:47] + "..." if ' ' in alignment.title else "No desc"
                
                print(f"{pdb_id:<10}{desc:<50}{hsp.expect:<12.2e}{identity_percentage:<15.2f}%", flush=True)
                
                if best_pdb_id is None:
                    best_pdb_id = pdb_id
                    best_similarity = identity_percentage
                    best_evalue = hsp.expect

    if best_pdb_id:
        print("\n" + "="*50, flush=True)
        print(f"📥 Fetching coordinate file for PDB ID: {best_pdb_id}...", flush=True)
        print("="*50, flush=True)
        
        pdbl = PDBList()
        downloaded_file = pdbl.retrieve_pdb_file(best_pdb_id, file_format="mmCif", pdir=".")
        print(f"\n🎉 File successfully saved to your folder as: {downloaded_file}", flush=True)

if __name__ == "__main__":
    # Test sequence (Human Insulin)
    raw_sequence = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    search_and_download_pdb(raw_sequence)
