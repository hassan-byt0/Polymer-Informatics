// Nodes
CREATE (:Polymer {
    id: $id,
    name: $name,
    psmiles: $psmiles,
    source: $source
})

CREATE (:Property {
    name: $name,
    value: $value,
    units: $units,
    conditions: $conditions
})

CREATE (:Literature {
    doi: $doi,
    title: $title,
    authors: $authors
})

// Relationships
MATCH (p:Polymer), (prop:Property)
CREATE (p)-[:HAS_PROPERTY {
    measurement_method: $method,
    confidence: $confidence
}]->(prop)

MATCH (p:Polymer), (lit:Literature)
CREATE (p)-[:CITED_IN]->(lit)
