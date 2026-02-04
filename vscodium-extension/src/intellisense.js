/**
 * IntelliSense Provider for Scientific APIs
 * Provides auto-completion, hover information, and documentation
 */

const vscode = require('vscode');

class ScientificAPICompletionProvider {
    constructor() {
        this.apis = this.loadScientificAPIs();
    }

    loadScientificAPIs() {
        return {
            'pubchempy': {
                module: 'pubchempy',
                alias: 'pcp',
                functions: [
                    {
                        name: 'get_compounds',
                        signature: 'get_compounds(identifier, namespace="cid", **kwargs)',
                        description: 'Retrieve Compound records from PubChem',
                        params: [
                            { name: 'identifier', type: 'str', description: 'Compound identifier' },
                            { name: 'namespace', type: 'str', description: 'Type of identifier: cid, name, smiles, inchi, inchikey, formula' }
                        ],
                        returns: 'List[Compound]',
                        example: 'compounds = pcp.get_compounds("aspirin", "name")'
                    },
                    {
                        name: 'get_properties',
                        signature: 'get_properties(properties, identifier, namespace="cid", **kwargs)',
                        description: 'Get specific properties for compounds',
                        params: [
                            { name: 'properties', type: 'List[str]', description: 'Properties to retrieve' },
                            { name: 'identifier', type: 'str', description: 'Compound identifier' }
                        ],
                        returns: 'List[Dict]',
                        example: 'props = pcp.get_properties(["MolecularFormula", "MolecularWeight"], "aspirin", "name")'
                    }
                ]
            },
            'biopython': {
                module: 'Bio.Entrez',
                alias: 'Entrez',
                functions: [
                    {
                        name: 'esearch',
                        signature: 'esearch(db, term, **kwargs)',
                        description: 'Search NCBI databases',
                        params: [
                            { name: 'db', type: 'str', description: 'Database to search: pubmed, gene, protein, nuccore' },
                            { name: 'term', type: 'str', description: 'Search query' },
                            { name: 'retmax', type: 'int', description: 'Maximum number of results', optional: true }
                        ],
                        returns: 'Handle',
                        example: 'handle = Entrez.esearch(db="pubmed", term="cancer", retmax=10)'
                    },
                    {
                        name: 'efetch',
                        signature: 'efetch(db, id, **kwargs)',
                        description: 'Fetch records from NCBI',
                        params: [
                            { name: 'db', type: 'str', description: 'Database name' },
                            { name: 'id', type: 'str', description: 'Record ID or list of IDs' },
                            { name: 'retmode', type: 'str', description: 'Return format: xml, text', optional: true }
                        ],
                        returns: 'Handle',
                        example: 'handle = Entrez.efetch(db="pubmed", id="12345", retmode="xml")'
                    }
                ]
            },
            'biothings_client': {
                module: 'biothings_client',
                alias: 'btc',
                functions: [
                    {
                        name: 'get_client',
                        signature: 'get_client(entity)',
                        description: 'Get a client for specific biological entity',
                        params: [
                            { name: 'entity', type: 'str', description: 'Entity type: gene, variant, chem, disease, taxon' }
                        ],
                        returns: 'Client',
                        example: 'mg = biothings_client.get_client("gene")'
                    }
                ]
            },
            'skbio': {
                module: 'skbio',
                alias: 'skbio',
                classes: [
                    {
                        name: 'DNA',
                        description: 'DNA sequence object',
                        methods: [
                            { name: 'complement', returns: 'DNA', description: 'Return complement sequence' },
                            { name: 'reverse_complement', returns: 'DNA', description: 'Return reverse complement' },
                            { name: 'gc_content', returns: 'float', description: 'Calculate GC content' },
                            { name: 'translate', returns: 'Protein', description: 'Translate to protein sequence' }
                        ],
                        example: 'seq = DNA("ACGTACGT")'
                    },
                    {
                        name: 'Protein',
                        description: 'Protein sequence object',
                        methods: [
                            { name: 'molecular_weight', returns: 'float', description: 'Calculate molecular weight' }
                        ],
                        example: 'prot = Protein("MKTAYIAKQRQ")'
                    }
                ]
            }
        };
    }

    provideCompletionItems(document, position, token, context) {
        const line = document.lineAt(position).text;
        const linePrefix = line.substring(0, position.character);

        const completions = [];

        // Detect import statements
        if (linePrefix.includes('import')) {
            for (const [key, api] of Object.entries(this.apis)) {
                const item = new vscode.CompletionItem(api.module, vscode.CompletionItemKind.Module);
                item.detail = `Scientific API: ${key}`;
                item.documentation = new vscode.MarkdownString(
                    `Import the ${key} scientific API module.\n\n` +
                    `Example: \`import ${api.module}\``
                );
                item.insertText = api.module;
                completions.push(item);
            }
        }

        // Detect function calls for specific modules
        for (const [key, api] of Object.entries(this.apis)) {
            if (api.functions && linePrefix.includes(api.alias + '.')) {
                api.functions.forEach(func => {
                    const item = new vscode.CompletionItem(func.name, vscode.CompletionItemKind.Function);
                    item.detail = func.signature;
                    item.documentation = new vscode.MarkdownString(
                        `**${func.description}**\n\n` +
                        `**Returns:** ${func.returns}\n\n` +
                        `**Example:**\n\`\`\`python\n${func.example}\n\`\`\``
                    );
                    
                    // Create snippet with parameters
                    const params = func.params.map((p, i) => `\${${i + 1}:${p.name}}`).join(', ');
                    item.insertText = new vscode.SnippetString(`${func.name}(${params})`);
                    
                    completions.push(item);
                });
            }

            if (api.classes && linePrefix.includes(api.alias + '.')) {
                api.classes.forEach(cls => {
                    const item = new vscode.CompletionItem(cls.name, vscode.CompletionItemKind.Class);
                    item.detail = `Class: ${cls.name}`;
                    item.documentation = new vscode.MarkdownString(
                        `**${cls.description}**\n\n` +
                        `**Example:**\n\`\`\`python\n${cls.example}\n\`\`\``
                    );
                    completions.push(item);
                });
            }
        }

        return completions;
    }

    provideHover(document, position, token) {
        const range = document.getWordRangeAtPosition(position);
        const word = document.getText(range);

        // Search for word in APIs
        for (const api of Object.values(this.apis)) {
            if (api.functions) {
                const func = api.functions.find(f => f.name === word);
                if (func) {
                    const markdown = new vscode.MarkdownString();
                    markdown.appendCodeblock(func.signature, 'python');
                    markdown.appendMarkdown(`\n${func.description}\n\n`);
                    markdown.appendMarkdown(`**Returns:** ${func.returns}\n\n`);
                    markdown.appendMarkdown(`**Example:**\n`);
                    markdown.appendCodeblock(func.example, 'python');
                    return new vscode.Hover(markdown);
                }
            }
        }

        return null;
    }

    provideSignatureHelp(document, position, token, context) {
        const line = document.lineAt(position).text;
        const beforeCursor = line.substring(0, position.character);
        
        // Find function name
        const match = beforeCursor.match(/(\w+)\s*\(/);
        if (!match) return null;

        const functionName = match[1];

        // Search for function in APIs
        for (const api of Object.values(this.apis)) {
            if (api.functions) {
                const func = api.functions.find(f => f.name === functionName);
                if (func) {
                    const sigHelp = new vscode.SignatureHelp();
                    const sig = new vscode.SignatureInformation(func.signature, func.description);
                    
                    func.params.forEach(param => {
                        const paramLabel = param.optional 
                            ? `${param.name}=${param.type}` 
                            : `${param.name}: ${param.type}`;
                        sig.parameters.push(new vscode.ParameterInformation(paramLabel, param.description));
                    });

                    sigHelp.signatures = [sig];
                    sigHelp.activeSignature = 0;
                    sigHelp.activeParameter = 0; // Could be calculated based on comma count
                    
                    return sigHelp;
                }
            }
        }

        return null;
    }
}

module.exports = ScientificAPICompletionProvider;
