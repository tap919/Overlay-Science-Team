#!/usr/bin/env python3
"""
Enhanced Scientific APIs Integration Module
Provides unified interface to multiple scientific databases and tools
"""

from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum

# API Integration Categories
class APICategory(Enum):
    COMPOUNDS = "compounds"
    LITERATURE = "literature"
    GENES = "genes"
    PROTEINS = "proteins"
    CANCER = "cancer"
    CHEMISTRY = "chemistry"
    BIOINFORMATICS = "bioinformatics"
    GENOMICS = "genomics"


@dataclass
class ScientificAPI:
    """Represents a scientific API with metadata"""
    name: str
    category: APICategory
    description: str
    install_command: str
    example_code: str
    documentation_url: str
    github_url: str
    requires_api_key: bool = False


class EnhancedScientificAPIs:
    """
    Enhanced scientific API integrations for VSCodium IDE
    Provides easy access to 15+ scientific databases and tools
    """
    
    def __init__(self):
        self.apis = self._initialize_apis()
    
    def _initialize_apis(self) -> Dict[str, ScientificAPI]:
        """Initialize all available scientific APIs"""
        return {
            # Compound & Drug APIs
            "pubchempy": ScientificAPI(
                name="PubChemPy",
                category=APICategory.COMPOUNDS,
                description="Wrapper for PubChem PUG REST API; fetches compounds by name/SMILES",
                install_command="pip install pubchempy",
                example_code="""import pubchempy as pcp
compounds = pcp.get_compounds('aspirin', 'name')
for c in compounds:
    print(f"{c.iupac_name}: {c.molecular_formula}")""",
                documentation_url="https://pubchempy.readthedocs.io",
                github_url="https://github.com/mcs07/PubChemPy"
            ),
            
            "chembl": ScientificAPI(
                name="ChEMBL Downloader",
                category=APICategory.CHEMISTRY,
                description="Downloads ChEMBL drug data for target interactions",
                install_command="pip install chembl-downloader",
                example_code="""from chembl_downloader import download
download('molecule', './chembl_data')""",
                documentation_url="https://chembl-downloader.readthedocs.io",
                github_url="https://github.com/chembl/chembl-downloader"
            ),
            
            "rdkit": ScientificAPI(
                name="RDKit",
                category=APICategory.CHEMISTRY,
                description="Cheminformatics and machine learning toolkit",
                install_command="pip install rdkit",
                example_code="""from rdkit import Chem
mol = Chem.MolFromSmiles('CCO')
print(Chem.MolToInchi(mol))""",
                documentation_url="https://www.rdkit.org/docs/",
                github_url="https://github.com/rdkit/rdkit"
            ),
            
            # Literature & Database APIs
            "biopython": ScientificAPI(
                name="Biopython",
                category=APICategory.LITERATURE,
                description="Entrez module for NCBI APIs (PubMed, Gene, Protein)",
                install_command="pip install biopython",
                example_code="""from Bio import Entrez
Entrez.email = 'your@email.com'
handle = Entrez.esearch(db='pubmed', term='cancer', retmax=10)
record = Entrez.read(handle)
print(record['IdList'])""",
                documentation_url="https://biopython.org",
                github_url="https://github.com/biopython/biopython"
            ),
            
            "metapub": ScientificAPI(
                name="MetaPub",
                category=APICategory.LITERATURE,
                description="Python wrapper for PubMed/NCBI with caching",
                install_command="pip install metapub",
                example_code="""from metapub import PubMedFetcher
fetch = PubMedFetcher()
article = fetch.article_by_pmid('12345')
print(article.title)""",
                documentation_url="https://metapub.readthedocs.io",
                github_url="https://github.com/nsh87/metapub"
            ),
            
            # Gene & Variant APIs
            "biothings": ScientificAPI(
                name="BioThings Client",
                category=APICategory.GENES,
                description="Client for BioThings APIs (MyGene.info, MyVariant.info)",
                install_command="pip install biothings-client",
                example_code="""from biothings_client import get_client
mg = get_client('gene')
results = mg.querymany(['TP53', 'BRCA1'], scopes='symbol')
for r in results:
    print(f"{r['symbol']}: {r['name']}")""",
                documentation_url="https://docs.biothings.io",
                github_url="https://github.com/biothings/biothings_client"
            ),
            
            "mygene": ScientificAPI(
                name="MyGene.info",
                category=APICategory.GENES,
                description="Gene annotation query service",
                install_command="pip install mygene",
                example_code="""import mygene
mg = mygene.MyGeneInfo()
gene = mg.getgene('1017', fields='symbol,name,pathway')
print(gene)""",
                documentation_url="https://docs.mygene.info",
                github_url="https://github.com/biothings/mygene.py"
            ),
            
            # Cancer Genomics APIs
            "pybioportal": ScientificAPI(
                name="PyBioPortal",
                category=APICategory.CANCER,
                description="Access to cBioPortal for cancer genomics data",
                install_command="pip install pybioportal",
                example_code="""from pybioportal import DataHub
dh = DataHub()
studies = dh.get_all_studies()
print(f"Available studies: {len(studies)}")""",
                documentation_url="https://pybioportal.readthedocs.io",
                github_url="https://github.com/asncd/pybioportal"
            ),
            
            # Protein & Structure APIs
            "biopandas": ScientificAPI(
                name="BioPandas",
                category=APICategory.PROTEINS,
                description="Working with molecular structures in pandas DataFrames",
                install_command="pip install biopandas",
                example_code="""from biopandas.pdb import PandasPdb
ppdb = PandasPdb().fetch_pdb('3eiy')
print(ppdb.df['ATOM'].head())""",
                documentation_url="https://biopandas.github.io/biopandas/",
                github_url="https://github.com/rasbt/biopandas"
            ),
            
            "prody": ScientificAPI(
                name="ProDy",
                category=APICategory.PROTEINS,
                description="Protein structural dynamics analysis",
                install_command="pip install prody",
                example_code="""import prody
protein = prody.parsePDB('1ubq')
print(protein)""",
                documentation_url="http://prody.csb.pitt.edu/",
                github_url="https://github.com/prody/ProDy"
            ),
            
            # Bioinformatics & Sequence APIs
            "skbio": ScientificAPI(
                name="scikit-bio",
                category=APICategory.BIOINFORMATICS,
                description="Bioinformatics tools for sequences, phylogenetics",
                install_command="pip install scikit-bio",
                example_code="""from skbio import DNA
seq = DNA('ACGTACGT')
print(f"GC Content: {seq.gc_content():.2%}")
print(f"Complement: {seq.complement()}")""",
                documentation_url="https://scikit.bio",
                github_url="https://github.com/biocore/scikit-bio"
            ),
            
            "pysam": ScientificAPI(
                name="pysam",
                category=APICategory.GENOMICS,
                description="Python wrapper for SAM/BAM/VCF/BCF files",
                install_command="pip install pysam",
                example_code="""import pysam
samfile = pysam.AlignmentFile('example.bam', 'rb')
for read in samfile.fetch():
    print(read)""",
                documentation_url="https://pysam.readthedocs.io",
                github_url="https://github.com/pysam-developers/pysam"
            ),
            
            # Genomics & Annotation APIs
            "biomart": ScientificAPI(
                name="BioMart",
                category=APICategory.GENOMICS,
                description="Wrapper for Ensembl BioMart API",
                install_command="pip install biomart",
                example_code="""from biomart import BiomartServer
server = BiomartServer('http://www.ensembl.org/biomart')
mart = server.datasets['hsapiens_gene_ensembl']
results = mart.search({'attributes': ['ensembl_gene_id', 'external_gene_name']})
for line in results.iter_lines():
    print(line.decode('utf-8'))""",
                documentation_url="https://pypi.org/project/biomart/",
                github_url="https://github.com/sebriois/biomart"
            ),
            
            "pyensembl": ScientificAPI(
                name="PyEnsembl",
                category=APICategory.GENOMICS,
                description="Python interface to Ensembl reference genome",
                install_command="pip install pyensembl",
                example_code="""import pyensembl
data = pyensembl.EnsemblRelease(75)
gene = data.gene_by_name('TP53')
print(f"Location: {gene.contig}:{gene.start}-{gene.end}")""",
                documentation_url="https://github.com/openvax/pyensembl",
                github_url="https://github.com/openvax/pyensembl"
            ),
            
            # Systems Biology & Pathway APIs
            "bioservices": ScientificAPI(
                name="BioServices",
                category=APICategory.BIOINFORMATICS,
                description="Access to 25+ biological web services",
                install_command="pip install bioservices",
                example_code="""from bioservices import KEGG
k = KEGG()
pathways = k.list('pathway', 'hsa')
print(pathways[:5])""",
                documentation_url="https://bioservices.readthedocs.io",
                github_url="https://github.com/cokelaer/bioservices"
            )
        }
    
    def get_by_category(self, category: APICategory) -> List[ScientificAPI]:
        """Get all APIs in a specific category"""
        return [api for api in self.apis.values() if api.category == category]
    
    def get_api(self, name: str) -> Optional[ScientificAPI]:
        """Get a specific API by name"""
        return self.apis.get(name.lower())
    
    def list_all_apis(self) -> List[str]:
        """List names of all available APIs"""
        return list(self.apis.keys())
    
    def generate_installation_script(self, api_names: List[str] = None) -> str:
        """Generate pip installation script for selected APIs"""
        if api_names is None:
            api_names = self.list_all_apis()
        
        commands = []
        for name in api_names:
            api = self.get_api(name)
            if api:
                commands.append(api.install_command)
        
        return "\n".join(commands)
    
    def get_example_notebook(self, category: APICategory = None) -> str:
        """Generate example Jupyter notebook code for a category"""
        if category:
            apis = self.get_by_category(category)
        else:
            apis = list(self.apis.values())[:5]
        
        notebook_code = "# Scientific API Examples\n\n"
        for api in apis:
            notebook_code += f"## {api.name}\n"
            notebook_code += f"# {api.description}\n"
            notebook_code += f"# Install: {api.install_command}\n\n"
            notebook_code += api.example_code + "\n\n"
        
        return notebook_code


# Convenience function for VSCodium extension
def get_enhanced_apis() -> EnhancedScientificAPIs:
    """Get instance of enhanced scientific APIs"""
    return EnhancedScientificAPIs()


if __name__ == "__main__":
    # Demo usage focused on automation bundles
    apis = EnhancedScientificAPIs()
    print("=== Enhanced Scientific APIs ===\n")
    print(f"Total APIs available: {len(apis.list_all_apis())}\n")
    for category in APICategory:
        category_apis = apis.get_by_category(category)
        print(f"{category.value.upper()}: {len(category_apis)} APIs")
    print("\n=== Open-Source Automation Bundles ===")
    bundles = {
        "literature_ingestion": {
            "purpose": "Fetch and normalize papers plus metadata automatically",
            "tools": "biopython, metapub, bioservices",
            "install": "pip install biopython metapub bioservices"
        },
        "compound_and_target_mining": {
            "purpose": "Rapidly search compounds, targets, and interactions for new hypotheses",
            "tools": "pubchempy, chembl-downloader, rdkit-pypi, biothings-client",
            "install": "pip install pubchempy chembl-downloader rdkit-pypi biothings-client"
        },
        "genomics_and_proteomics": {
            "purpose": "Analyze sequences, variants, and protein structures in automated pipelines",
            "tools": "pysam, scikit-bio, pybioportal, biopandas, prody",
            "install": "pip install pysam scikit-bio pybioportal biopandas prody"
        },
        "analytics_and_insights": {
            "purpose": "Turn raw outputs into dashboards and insight reports automatically",
            "tools": "pandas, scipy, scikit-learn, plotly",
            "install": "pip install pandas scipy scikit-learn plotly"
        }
    }
    for name, bundle in bundles.items():
        print(f"\n[{name}]")
        print(f"Purpose : {bundle['purpose']}")
        print(f"Tools   : {bundle['tools']}")
        print(f"Install : {bundle['install']}")
    print("\n=== Quickstart Script ===")
    print("#!/usr/bin/env bash")
    print("set -eo pipefail")
    print("echo \"Installing open-source science automation stack...\"")
    for bundle in bundles.values():
        print(f"echo \"- {bundle['purpose']}\"")
        print(bundle["install"])
    print("echo \"Creating starter notebooks...\"")
    print("python examples/protein-folding-study/study_setup.py 2>/dev/null || true")
    print("echo \"Done. Launch your science lab!\"")
