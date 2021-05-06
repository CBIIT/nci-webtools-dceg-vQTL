import React, { useContext, useState, useEffect, useRef } from 'react';
import { RootContext } from '../../../index';
import {
  Form,
  Button,
  Row,
  Col,
  Popover,
  OverlayTrigger,
} from 'react-bootstrap';
import Overlay from 'react-bootstrap/Overlay';
import Tooltip from 'react-bootstrap/Tooltip';
import { useDispatch, useSelector } from 'react-redux';
import {
  qtlsGWASCalculation,
  uploadFile,
  updateQTLsGWAS,
  getPublicGTEx,
  updateAlert,
  submitQueue,
} from '../../../services/actions';
import Select from '../../controls/select/select';
import Accordions from '../../controls/accordions/accordions';
import { PopulationSelect } from '../../controls/population-select/population-select';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faInfoCircle } from '@fortawesome/free-solid-svg-icons';
const { v1: uuidv1 } = require('uuid');

export function QTLsGWASForm() {
  const dispatch = useDispatch();
  const { getInitialState } = useContext(RootContext);

  const quantificationFileControl = useRef(null);
  const genotypeFileControl = useRef(null);

  const [_associationFile, _setAssociationFile] = useState('');
  const [_quantificationFile, _setQuantificationFile] = useState('');
  const [_genotypeFile, _setGenotypeFile] = useState('');
  const [_LDFile, _setLDFile] = useState('');
  const [_gwasFile, _setGwasFile] = useState('');
  const [showQuantificationTooltip, setShowQuantificationTooltip] = useState(
    false
  );
  const [showGenotypeTooltip, setShowGenotypeTooltip] = useState(false);
  const [tissueOnly, viewTissueOnly] = useState(false);
  const [phenotypeOnly, viewPhenotypeOnly] = useState(false);
  const [showAdditionalInput, setAdditionalInput] = useState(false);

  const {
    select_qtls_samples,
    select_gwas_sample,
    associationFile,
    quantificationFile,
    genotypeFile,
    gwasFile,
    LDFile,
    request,
    select_pop,
    select_gene,
    select_dist,
    select_ref,
    recalculateAttempt,
    recalculatePop,
    recalculateGene,
    recalculateDist,
    recalculateRef,
    submitted,
    isLoading,
    isError,
    publicGTEx,
    publicLoading,
    genome,
    genomeOptions,
    qtlProject,
    qtlProjectOptions,
    gwasProject,
    gwasProjectOptions,
    ldProject,
    ldProjectOptions,
    xQtl,
    xQtlOptions,
    tissue,
    tissueOptions,
    phenotype,
    phenotypeOptions,
    select_chromosome,
    select_position,
    isQueue,
    email,
    qtlPublic,
    gwasPublic,
    ldPublic,
    qtlKey,
    ldKey,
    gwasKey,
  } = useSelector((state) => state.qtlsGWAS);

  useEffect(() => _setAssociationFile(associationFile), [associationFile]);
  useEffect(() => _setQuantificationFile(quantificationFile), [
    quantificationFile,
  ]);
  useEffect(() => _setGenotypeFile(genotypeFile), [genotypeFile]);
  useEffect(() => _setLDFile(LDFile), [LDFile]);
  useEffect(() => _setGwasFile(gwasFile), [gwasFile]);
  useEffect(() => {
    if (select_qtls_samples || select_gwas_sample) setAdditionalInput(true);
    else setAdditionalInput(false);
  }, [select_qtls_samples, select_gwas_sample]);
  useEffect(() => {
    if (!Object.keys(publicGTEx).length) dispatch(getPublicGTEx());
  }, [publicGTEx]);
  // inital population of public params
  useEffect(() => {
    if (Object.keys(publicGTEx).length && !genomeOptions.length)
      getGenomeOptions();
  }, [publicGTEx, genomeOptions]);
  useEffect(() => {
    if (Object.keys(publicGTEx).length && genome.value)
      populatePublicParameters();
  }, [genome.value]);
  useEffect(() => {
    if (Object.keys(publicGTEx).length && qtlProject && qtlPublic)
      handleQtlProject(qtlProject);
  }, [qtlPublic, tissueOnly]);
  useEffect(() => {
    if (Object.keys(publicGTEx).length && gwasProject && gwasPublic)
      handleGwasProject(gwasProject);
  }, [gwasPublic, phenotypeOnly]);
  useEffect(() => {
    if (Object.keys(publicGTEx).length && ldProject && ldPublic)
      handleLdProject(ldProject);
  }, [ldPublic]);

  function getGenomeOptions() {
    const data = publicGTEx['cis-QTL dataset'];
    const genomeOptions = [
      ...new Set(data.map((row) => row.Genome_build)),
    ].map((genome) => ({ value: genome, label: genome }));

    dispatch(
      updateQTLsGWAS({
        genome: genomeOptions[0],
        genomeOptions: genomeOptions,
      })
    );
  }

  function getProjectOptions(data) {
    return [
      ...new Set(
        data
          .filter((row) => row['Genome_build'] == genome.value)
          .map((row) => row['Project'])
      ),
    ].map((project) => ({ value: project, label: project }));
  }

  function getXqtlOptions(data, project) {
    return [
      ...new Set(
        data
          .filter(
            (row) => row.Genome_build == genome.value && row.Project == project
          )
          .map((row) => row.xQTL)
      ),
    ].map((xQtl) => ({ value: xQtl, label: xQtl }));
  }

  function getTissueOptions(data, project, xQtl) {
    return !tissueOnly
      ? data
          .filter(
            (row) =>
              row['Genome_build'] == genome.value &&
              row['Project'] == project &&
              row['xQTL'] == xQtl
          )
          .map((row) => ({ value: row.Tissue, label: row.Tissue }))
      : data
          .filter((row) => row['Genome_build'] == genome.value)
          .map((row) => ({ value: row.Tissue, label: row.Full_Name }));
  }

  function getPhenotypeOptions(data, project) {
    return !phenotypeOnly
      ? data
          .filter(
            (row) =>
              row['Genome_build'] == genome.value && row['Project'] == project
          )
          .map((row) => ({ value: row.Phenotype, label: row.Phenotype }))
      : data
          .filter((row) => row['Genome_build'] == genome.value)
          .map((row) => ({ value: row.Phenotype, label: row.Full_Name }));
  }

  function populatePublicParameters() {
    const qtlData = publicGTEx['cis-QTL dataset'];
    const ldData = publicGTEx['LD dataset'];
    const gwasData = publicGTEx['GWAS dataset'];

    const qtlProjectOptions = getProjectOptions(qtlData);
    const xQtlOptions = getXqtlOptions(qtlData, qtlProjectOptions[0].value);
    const tissueOptions = getTissueOptions(
      qtlData,
      qtlProjectOptions[0].value,
      xQtlOptions[0].value
    );

    const ldProjectOptions = getProjectOptions(ldData);

    const gwasProjectOptions = getProjectOptions(gwasData);
    const phenotypeOptions = gwasProjectOptions.length
      ? getPhenotypeOptions(gwasData, gwasProjectOptions[0].value)
      : [];

    dispatch(
      updateQTLsGWAS({
        qtlProject: qtlProjectOptions[0],
        ldProject: ldProjectOptions[0],
        xQtl: xQtlOptions[0],
        tissue: tissueOptions[0],
        gwasProject: gwasProjectOptions[0],
        phenotype: phenotypeOptions[0] || '',

        qtlProjectOptions: qtlProjectOptions,
        ldProjectOptions: ldProjectOptions,
        xQtlOptions: xQtlOptions,
        tissueOptions: tissueOptions,
        gwasProjectOptions: gwasProjectOptions,
        phenotypeOptions: phenotypeOptions || [],
      })
    );
  }

  function handleQtlProject(project) {
    const data = publicGTEx['cis-QTL dataset'];
    const xQtlOptions = getXqtlOptions(data, project.value);
    const tissueOptions = getTissueOptions(
      data,
      project.value,
      xQtlOptions[0].value
    );
    const qtlKey = data
      .filter(
        (row) =>
          row.Genome_build == genome.value &&
          row.Project == project.value &&
          row.xQTL == xQtlOptions[0].value &&
          row.Tissue == tissueOptions[0].value
      )[0]
      .Biowulf_full_path.replace('/data/Brown_lab/ZTW_KB_Datasets/vQTL2/', '');

    dispatch(
      updateQTLsGWAS({
        qtlProject: project,
        xQtl: xQtlOptions[0],
        xQtlOptions: xQtlOptions,
        tissue: tissueOptions[0],
        tissueOptions: tissueOptions,
        qtlKey: qtlKey,
      })
    );
  }

  function handleGwasProject(project) {
    const data = publicGTEx['GWAS dataset'];
    const phenotypeOptions = getPhenotypeOptions(data, project.value);
    const gwasKey = data
      .filter(
        (row) =>
          row.Genome_build == genome.value &&
          row.Project == project.value &&
          row.Phenotype == phenotypeOptions[0].value
      )[0]
      .Biowulf_full_path.replace('/data/Brown_lab/ZTW_KB_Datasets/vQTL2/', '');

    dispatch(
      updateQTLsGWAS({
        gwasProject: project,
        phenotype: phenotypeOptions[0],
        phenotypeOptions: phenotypeOptions,
        gwasKey: gwasKey,
      })
    );
  }

  function handleLdProject(project) {
    const ldKey =
      project.value == '1000genome'
        ? publicGTEx['LD dataset']
            .filter(
              (row) =>
                row.Genome_build == genome.value &&
                row.Project == project.value &&
                row.Chromosome == select_chromosome.value
            )[0]
            .Biowulf_full_path.replace(
              '/data/Brown_lab/ZTW_KB_Datasets/vQTL2/',
              ''
            )
        : true;

    dispatch(updateQTLsGWAS({ ldProject: project, ldKey: ldKey }));
  }

  // qtl type
  function handleXqtl(xQtl) {
    const data = publicGTEx['cis-QTL dataset'];
    const tissueOptions = getTissueOptions(data, qtlProject.value, xQtl.value);
    const qtlKey = data
      .filter(
        (row) =>
          row.Genome_build == genome.value &&
          row.Project == qtlProject.value &&
          row.xQTL == xQtl.value &&
          row.Tissue == tissueOptions[0].value
      )[0]
      .Biowulf_full_path.replace('/data/Brown_lab/ZTW_KB_Datasets/vQTL2/', '');

    dispatch(
      updateQTLsGWAS({
        xQtl: xQtl,
        tissue: tissueOptions[0],
        tissueOptions: tissueOptions,
        qtlKey: qtlKey,
      })
    );
  }

  function handleTissue(tissue) {
    const qtlKey = publicGTEx['cis-QTL dataset']
      .filter(
        (row) =>
          row.Genome_build == genome.value &&
          row.Project == qtlProject.value &&
          row.xQTL == xQtl.value &&
          row.Tissue == tissue.value
      )[0]
      .Biowulf_full_path.replace('/data/Brown_lab/ZTW_KB_Datasets/vQTL2/', '');
    dispatch(updateQTLsGWAS({ tissue: tissue, qtlKey: qtlKey }));
  }

  function handlePhenotype(phenotype) {
    const gwasKey = publicGTEx['GWAS dataset']
      .filter(
        (row) =>
          row.Genome_build == genome.value &&
          row.Project == gwasProject.value &&
          row.Phenotype == phenotype.value
      )[0]
      .Biowulf_full_path.replace('/data/Brown_lab/ZTW_KB_Datasets/vQTL2/', '');

    dispatch(
      updateQTLsGWAS({
        phenotype: phenotype,
        gwasKey: gwasKey,
      })
    );
  }

  const handleReset = () => {
    window.location.hash = '#/qtls';

    dispatch(
      updateQTLsGWAS({ ...getInitialState().qtlsGWAS, publicGTEx: publicGTEx })
    );
    dispatch(updateAlert(getInitialState().alert));
    _setAssociationFile('');
    _setQuantificationFile('');
    _setGenotypeFile('');
    _setLDFile('');
    _setGwasFile('');
    viewTissueOnly(false);
    viewPhenotypeOnly(false);
  };

  async function handleSubmit() {
    if (_quantificationFile && !_genotypeFile) {
      setShowGenotypeTooltip(true);
      setTimeout(() => setShowGenotypeTooltip(false), 5000);
      return;
    }
    if (!_quantificationFile && _genotypeFile) {
      setShowQuantificationTooltip(true);
      setTimeout(() => setShowQuantificationTooltip(false), 5000);
      return;
    }
    const request = uuidv1();

    await dispatch(
      uploadFile({
        // dataFiles: [
        //   _associationFile,
        //   _quantificationFile,
        //   _genotypeFile,
        //   _LDFile,
        //   _gwasFile,
        // ],
        associationFile: _associationFile,
        quantificationFile: _quantificationFile,
        genotypeFile: _genotypeFile,
        LDFile: _LDFile,
        gwasFile: _gwasFile,
        associationFileName: _associationFile ? _associationFile.name : false,
        quantificationFileName: _quantificationFile
          ? _quantificationFile.name
          : false,
        genotypeFileName: _genotypeFile ? _genotypeFile.name : false,
        LDFileName: _LDFile ? _LDFile.name : false,
        gwasFileName: _gwasFile ? _gwasFile.name : false,
        request,
      })
    );

    const params = {
      request,
      select_qtls_samples,
      select_gwas_sample,
      associationFile: (_associationFile && _associationFile.name) || false,
      quantificationFile:
        (_quantificationFile && _quantificationFile.name) || false,
      genotypeFile: (_genotypeFile && _genotypeFile.name) || false,
      gwasFile: (_gwasFile && _gwasFile.name) || false,
      LDFile: (_LDFile && _LDFile.name) || false,
      select_pop,
      select_gene,
      select_dist,
      select_ref,
      recalculateAttempt,
      recalculatePop,
      recalculateGene,
      recalculateDist,
      recalculateRef,
      ldProject: ldProject.value,
      qtlKey: qtlPublic ? qtlKey : false,
      ldKey: ldPublic ? ldKey : false,
      gwasKey: gwasPublic ? gwasKey : false,
      select_chromosome: select_chromosome.value,
      select_position: select_position,
      email: email,
    };

    if (isQueue) {
      dispatch(submitQueue(params));
    } else {
      dispatch(qtlsGWASCalculation(params));
    }
  }

  const accordionComponents = [
    {
      title: 'QTLs Data',
      component: (
        <>
          <Row>
            <Form.Group className="col-sm-12">
              <div className="d-flex">
                <Form.Label className="mb-0 mr-auto">
                  Association (QTL) Data <span style={{ color: 'red' }}>*</span>
                </Form.Label>
                <Form.Check
                  title="Association (QTL) Public Data Checkbox"
                  disabled={submitted || select_qtls_samples}
                  inline
                  id="qtlSource"
                  label="Public"
                  type="checkbox"
                  checked={qtlPublic}
                  onChange={(_) => {
                    dispatch(
                      updateQTLsGWAS({
                        qtlPublic: !qtlPublic,
                        ...(!qtlPublic && { select_ref: false }),
                      })
                    );
                    _setAssociationFile('');
                  }}
                />
              </div>
              {qtlPublic ? (
                <div className="mt-2">
                  <Row>
                    <Col>
                      <Form.Check
                        title="Association (QTL) Public Data Tissue Only Checkbox"
                        id="tissueOnly"
                        label="Tissue Only"
                        type="checkbox"
                        disabled={
                          submitted ||
                          publicLoading ||
                          !Object.keys(publicGTEx).length
                        }
                        checked={tissueOnly}
                        onChange={(_) => {
                          viewTissueOnly(!tissueOnly);
                        }}
                      />
                    </Col>
                  </Row>
                  {!tissueOnly && (
                    <>
                      <Row>
                        <Col>
                          <Select
                            disabled={publicLoading || tissueOnly || submitted}
                            id="project"
                            label="Project"
                            value={qtlProject}
                            options={qtlProjectOptions}
                            onChange={handleQtlProject}
                          />
                        </Col>
                      </Row>
                      <Row>
                        <Col>
                          <Select
                            disabled={publicLoading || tissueOnly || submitted}
                            id="qtlType"
                            label="QTL Type"
                            value={xQtl}
                            options={xQtlOptions}
                            onChange={handleXqtl}
                          />
                        </Col>
                      </Row>
                    </>
                  )}
                  <Row>
                    <Col>
                      <Select
                        disabled={publicLoading || submitted}
                        id="tissue"
                        label="Tissue"
                        value={tissue}
                        options={tissueOptions}
                        onChange={handleTissue}
                      />
                    </Col>
                  </Row>
                </div>
              ) : (
                <Form.File
                  title="Association (QTL) Data User File Upload Input"
                  id="qtls-association-file"
                  disabled={submitted || select_qtls_samples}
                  key={_associationFile}
                  label={
                    _associationFile
                      ? _associationFile.name ||
                        _associationFile.filename ||
                        _associationFile
                      : select_qtls_samples
                      ? 'MX2.eQTL.txt'
                      : 'Choose File'
                  }
                  onChange={(e) => {
                    _setAssociationFile(e.target.files[0]);
                  }}
                  // accept=".tsv, .txt"
                  // isInvalid={checkValid ? !validFile : false}
                  // feedback="Please upload a data file"
                  // onChange={(e) => {
                  //     setInput(e.target.files[0]);
                  //     mergeVisualize({
                  //     storeFilename: e.target.files[0].name,
                  //     });
                  // }}
                  custom
                />
              )}
            </Form.Group>
            <Form.Group className="col-sm-12">
              <div className="d-flex">
                <Form.Label className="mb-0 mr-auto">GWAS Data</Form.Label>
                <Form.Check
                  title="GWAS Public Data Checkbox"
                  disabled={submitted || select_qtls_samples}
                  inline
                  id="gwasSource"
                  label="Public"
                  type="checkbox"
                  checked={gwasPublic}
                  onChange={(_) => {
                    dispatch(
                      updateQTLsGWAS({
                        gwasPublic: !gwasPublic,
                        ...(!ldPublic && { select_ref: false }),
                      })
                    );
                    _setGwasFile('');
                  }}
                />
              </div>
              {gwasPublic ? (
                <div className="mt-2">
                  <Row>
                    <Col>
                      <Form.Check
                        title="GWAS Public Data Phenotype Only Checkbox"
                        id="phenotypeOnly"
                        label="Phenotype Only"
                        type="checkbox"
                        disabled={
                          submitted ||
                          publicLoading ||
                          !Object.keys(publicGTEx).length
                        }
                        checked={phenotypeOnly}
                        onChange={(_) => {
                          viewPhenotypeOnly(!phenotypeOnly);
                        }}
                      />
                    </Col>
                  </Row>
                  {!phenotypeOnly && (
                    <Row>
                      <Col>
                        <Select
                          disabled={
                            submitted ||
                            publicLoading ||
                            phenotypeOnly ||
                            !gwasProjectOptions.length
                          }
                          id="gwasProject"
                          label="Project"
                          value={gwasProject}
                          options={gwasProjectOptions}
                          onChange={handleGwasProject}
                        />
                      </Col>
                    </Row>
                  )}

                  <Row>
                    <Col>
                      <Select
                        disabled={
                          submitted || publicLoading || !phenotypeOptions.length
                        }
                        id="gwasPhenotype"
                        label="Phenotype"
                        value={phenotype}
                        options={phenotypeOptions}
                        onChange={handlePhenotype}
                      />
                    </Col>
                  </Row>
                </div>
              ) : (
                <Form.File
                  title="GWAS Data User File Upload Input"
                  id="qtls-gwas-file"
                  disabled={submitted || select_gwas_sample}
                  key={_gwasFile}
                  label={
                    _gwasFile
                      ? _gwasFile.name || _gwasFile.filename || _gwasFile
                      : select_gwas_sample
                      ? 'MX2.GWAS.rs.txt'
                      : 'Choose File'
                  }
                  onChange={(e) => {
                    _setGwasFile(e.target.files[0]);
                  }}
                  // accept=".tsv, .txt"
                  // isInvalid={checkValid ? !validFile : false}
                  // feedback="Please upload a data file"
                  // onChange={(e) => {
                  //     setInput(e.target.files[0]);
                  //     mergeVisualize({
                  //     storeFilename: e.target.files[0].name,
                  //     });
                  // }}
                  custom
                />
              )}
              <Form.Text className="text-muted">
                <i>Upload locus specific region, &le; 5Mb size</i>
              </Form.Text>
            </Form.Group>
          </Row>
        </>
      ),
    },
    {
      title: 'Locus Quantification',
      collapseDefault: true,
      component: (
        <>
          <Row>
            <Form.Group className="col-sm-12">
              <Form.Label className="mb-0">Quantification Data</Form.Label>
              <Form.File
                title="Quantification Data User File Upload Input"
                ref={quantificationFileControl}
                id="qtls-quantification-file"
                disabled={submitted || select_qtls_samples}
                key={_quantificationFile}
                label={
                  _quantificationFile
                    ? _quantificationFile.name ||
                      _quantificationFile.filename ||
                      _quantificationFile
                    : select_qtls_samples
                    ? 'MX2.quantification.txt'
                    : 'Choose File'
                }
                onChange={(e) => {
                  _setQuantificationFile(e.target.files[0]);
                }}
                // accept=".tsv, .txt"
                // isInvalid={checkValid ? !validFile : false}
                // feedback="Please upload a data file"
                // onChange={(e) => {
                //     setInput(e.target.files[0]);
                //     mergeVisualize({
                //     storeFilename: e.target.files[0].name,
                //     });
                // }}
                custom
              />
              <Overlay
                target={quantificationFileControl.current}
                show={showQuantificationTooltip}
                placement="bottom"
              >
                {(props) => (
                  <Tooltip id="overlay-example" {...props}>
                    Please input accompanying Quantification Data File with
                    Genotype Data File.
                  </Tooltip>
                )}
              </Overlay>
            </Form.Group>
            <Form.Group className="col-sm-12">
              <Form.Label className="mb-0">Genotype Data</Form.Label>
              <Form.File
                title="Genotype Data User File Upload Input"
                ref={genotypeFileControl}
                id="qtls-genotype-file"
                disabled={submitted || select_qtls_samples}
                key={_genotypeFile}
                label={
                  _genotypeFile
                    ? _genotypeFile.name ||
                      _genotypeFile.filename ||
                      _genotypeFile
                    : select_qtls_samples
                    ? 'MX2.genotyping.txt'
                    : 'Choose File'
                }
                onChange={(e) => {
                  _setGenotypeFile(e.target.files[0]);
                }}
                // accept=".tsv, .txt"
                // isInvalid={checkValid ? !validFile : false}
                // feedback="Please upload a data file"
                // onChange={(e) => {
                //     setInput(e.target.files[0]);
                //     mergeVisualize({
                //     storeFilename: e.target.files[0].name,
                //     });
                // }}
                custom
              />
              <Overlay
                target={genotypeFileControl.current}
                show={showGenotypeTooltip}
                placement="bottom"
              >
                {(props) => (
                  <Tooltip id="overlay-example" {...props}>
                    Please input accompanying Genotype Data File with
                    Quantification Data File.
                  </Tooltip>
                )}
              </Overlay>
            </Form.Group>
          </Row>
        </>
      ),
    },
    {
      title: 'LD Information',
      component: (
        <Row>
          <Form.Group className="col-sm-12 mb-0">
            <div className="d-flex">
              <Form.Label className="mb-0 mr-auto">
                LD Data{' '}
                <OverlayTrigger
                  trigger="click"
                  placement="top"
                  overlay={
                    <Popover id="popover-basic">
                      <Popover.Title as="h3">LD Information</Popover.Title>
                      <Popover.Content>
                        <p>Default: 1KG Phase 3, EUR</p>
                      </Popover.Content>
                    </Popover>
                  }
                  rootClose
                >
                  <Button
                    variant="link"
                    className="p-0 font-weight-bold"
                    aria-label="LD Information additional info"
                  >
                    <FontAwesomeIcon icon={faInfoCircle} />
                  </Button>
                </OverlayTrigger>
              </Form.Label>
              <Form.Check
                disabled={submitted || select_qtls_samples}
                inline
                id="ldSource"
                label="Public"
                type="checkbox"
                checked={ldPublic}
                onChange={(_) => {
                  dispatch(
                    updateQTLsGWAS({
                      ldPublic: !ldPublic,
                      select_pop: false,
                      ...(!ldPublic && { select_ref: false }),
                    })
                  );
                  _setLDFile('');
                }}
              />
            </div>
            {ldPublic ? (
              <div>
                <Form.Row>
                  <Col>
                    <Select
                      disabled={publicLoading || submitted}
                      id="ldProject"
                      label="Project"
                      value={ldProject}
                      options={ldProjectOptions}
                      onChange={handleLdProject}
                    />
                  </Col>
                </Form.Row>
                <Form.Row>
                  <Col>
                    <Select
                      disabled={submitted}
                      id="chromosome"
                      label="Chromosome"
                      value={select_chromosome}
                      options={[
                        ...Array.from({ length: 22 }, (_, i) => ({
                          value: i + 1,
                          label: i + 1,
                        })),
                        {
                          value: 'X',
                          label: 'X',
                        },
                        {
                          value: 'Y',
                          label: 'Y',
                        },
                      ]}
                      onChange={(chromosome) => {
                        dispatch(
                          updateQTLsGWAS({ select_chromosome: chromosome })
                        );
                      }}
                    />
                  </Col>
                </Form.Row>
                <Form.Row>
                  <Col>
                    <Form.Label className="mb-0">
                      Population <span style={{ color: 'red' }}>*</span>{' '}
                    </Form.Label>
                    <PopulationSelect
                      id="qtls-results-population-input"
                      disabled={submitted || !ldPublic}
                    />
                  </Col>
                </Form.Row>
              </div>
            ) : (
              <>
                <Form.File
                  id="qtls-ld-file"
                  disabled={submitted || select_qtls_samples}
                  key={_LDFile}
                  label={
                    _LDFile
                      ? _LDFile.name || _LDFile.filename || _LDFile
                      : select_qtls_samples
                      ? 'MX2.LD.gz'
                      : 'Choose File'
                  }
                  onChange={(e) => {
                    _setLDFile(e.target.files[0]);
                  }}
                  // accept=".tsv, .txt"
                  // isInvalid={checkValid ? !validFile : false}
                  // feedback="Please upload a data file"
                  // onChange={(e) => {
                  //     setInput(e.target.files[0]);
                  //     mergeVisualize({
                  //     storeFilename: e.target.files[0].name,
                  //     });
                  // }}
                  custom
                />
              </>
            )}
          </Form.Group>
        </Row>
      ),
    },
    {
      title: 'Locus Information',
      component: (
        <>
          {' '}
          <Row>
            <Col>
              <Form.Label className="mb-0">
                cis-QTL Distance <span style={{ color: 'red' }}>*</span>{' '}
                <small>
                  <i>(+/- Kb up to 5Mb)</i>
                </small>
              </Form.Label>
              <Form.Control
                title="cis-QTL Distance Input"
                aria-label="cis-QTL Distance Input"
                type="number"
                min="1"
                max="2000"
                id="qtls-distance-input"
                disabled={submitted}
                onChange={(e) => {
                  dispatch(updateQTLsGWAS({ select_dist: e.target.value }));
                }}
                value={select_dist}
                isInvalid={select_dist < 1 || select_dist > 200}
                // custom
              />
              <Form.Control.Feedback type="invalid">
                Enter distance between 1 and 200Kb.
              </Form.Control.Feedback>
            </Col>
          </Row>
          {qtlPublic || ldPublic || gwasPublic ? (
            <Row>
              <Col>
                <Form.Row>
                  <Col>
                    <Form.Label className="mb-0">Position</Form.Label>
                    <Form.Control
                      title="LD Reference Position Input"
                      aria-label="LD Reference Position Input"
                      id="select_position"
                      disabled={submitted}
                      onChange={(e) => {
                        dispatch(
                          updateQTLsGWAS({ select_position: e.target.value })
                        );
                      }}
                      placeholder="e.g. 100000"
                      value={select_position}
                    />
                  </Col>
                </Form.Row>
              </Col>
            </Row>
          ) : (
            <Row>
              <Col>
                <Form.Label className="mb-0">
                  SNP{' '}
                  <small>
                    <i>(Default: lowest GWAS P-value SNP)</i>
                  </small>
                </Form.Label>
                <Form.Control
                  title="LD Reference SNP Input"
                  aria-label="LD Reference SNP Input"
                  id="qtls-snp-input"
                  disabled={submitted}
                  onChange={(e) => {
                    dispatch(updateQTLsGWAS({ select_ref: e.target.value }));
                  }}
                  f
                  value={select_ref ? select_ref : ''}
                  isInvalid={
                    select_ref &&
                    select_ref.length > 0 &&
                    !/^rs\d+$/.test(select_ref)
                  }
                  // custom
                />
                <Form.Control.Feedback type="invalid">
                  Enter valid RS number. Leave empty for default.
                </Form.Control.Feedback>
              </Col>
            </Row>
          )}
        </>
      ),
    },
    {
      title: 'Queue Job',
      component: (
        <div>
          <Row>
            <Col>
              <Form.Group controlId="toggleQueue" className="mb-0">
                <Form.Check inline>
                  <Form.Check.Input
                    className="mr-0"
                    title="Choose Queued Submission Checkbox"
                    type="checkbox"
                    disabled={submitted}
                    checked={isQueue}
                    onChange={(_) => {
                      dispatch(updateQTLsGWAS({ isQueue: !isQueue }));
                    }}
                  />
                </Form.Check>
                <Form.Label className="mr-auto">
                  Submit this job to a Queue
                </Form.Label>
              </Form.Group>
            </Col>
          </Row>
          <Row>
            <Col>
              <Form.Group controlId="email">
                <Form.Control
                  title="Queued Submission Email Input"
                  aria-label="Queued Submission Email Input"
                  placeholder="Enter Email"
                  size="sm"
                  value={email}
                  type="email"
                  onChange={(e) =>
                    dispatch(updateQTLsGWAS({ email: e.target.value }))
                  }
                  disabled={!isQueue || submitted}
                  // isInvalid={isQueue && checkValid ? !validEmail : false}
                />
                <Form.Control.Feedback type="invalid">
                  Please provide a valid email
                </Form.Control.Feedback>
                <Form.Text className="text-muted">
                  <i>
                    Note: if sending to queue, when computation is completed, a
                    notification will be sent to the e-mail entered above.
                  </i>
                </Form.Text>
              </Form.Group>
            </Col>
          </Row>
        </div>
      ),
    },
  ];

  return (
    <Form>
      <Row>
        <Col>
          <Select
            className="border rounded p-2"
            disabled={!genomeOptions.length || submitted}
            id="genomeBuild"
            label="Genome Build"
            value={genome}
            options={genomeOptions}
            onChange={(genome) => dispatch(updateQTLsGWAS({ genome: genome }))}
          />
        </Col>
      </Row>
      <Row>
        <Col sm="6">
          <Button
            className="p-0 font-14"
            variant="link"
            onClick={() => {
              _setGwasFile('');
              if (select_qtls_samples) {
                dispatch(
                  updateQTLsGWAS({
                    select_qtls_samples: true,
                    select_gwas_sample: true,
                    qtlPublic: false,
                    gwasPublic: false,
                    ldPublic: false,
                    select_pop: false,
                  })
                );
              } else {
                dispatch(
                  updateQTLsGWAS({
                    select_qtls_samples: false,
                    select_gwas_sample: false,
                  })
                );
              }
            }}
            disabled={submitted}
          >
            {!select_gwas_sample ? 'Load' : 'Unload'} Sample Data
          </Button>
        </Col>

        <Col sm="6">
          <a className="font-14" href="assets/files/MX2.examples.gz" download>
            Download Sample Data
          </a>
        </Col>
      </Row>
      <Accordions components={accordionComponents} bodyClass="p-2" />
      <Row className="mt-2">
        <Col sm="6">
          <Button
            // disabled={loading.active}
            className="w-100"
            variant={isError ? 'danger' : 'secondary'}
            onClick={() => handleReset()}
            disabled={submitted && isLoading}
          >
            Reset
          </Button>
        </Col>
        <Col sm="6">
          <Button
            // disabled={submitted || loading.active}
            className="w-100"
            variant="primary"
            type="button"
            onClick={() => {
              handleSubmit();
            }}
            disabled={
              submitted ||
              (!_associationFile && !select_qtls_samples && !qtlPublic) ||
              select_dist.length <= 0 ||
              select_dist < 1 ||
              select_dist > 200 ||
              (select_ref &&
                select_ref.length > 0 &&
                !/^rs\d+$/.test(select_ref))
              //   ||
              // (ldPublic && (!select_pop || select_pop.length <= 0))
            }
          >
            {isQueue ? 'Submit' : 'Calculate'}
          </Button>
        </Col>
      </Row>
    </Form>
  );
}
