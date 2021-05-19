import React, { useState, useEffect } from 'react';
import { useDispatch, useSelector } from 'react-redux';
// import { actions } from '../../../../services/store';
import { updateQTLsGWAS, fetchResults } from '../../../../services/actions';
import { Tabs, Tab, Form, Button } from 'react-bootstrap';
import { LocusQC } from './locus-qc';
import { LocusLD } from './locus-ld';
import { LocusAlignment } from './locus-alignment';
import { LocusColocalization } from './locus-colocalization';
import { LocusQuantifiation } from './locus-quantification';
import { LocusTable } from './locus-table';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { LocusDownload } from './locus-download';
import { QTLsGWASResultsForm } from './qtls-gwas-results-form';

export function QTLsGWASResults({ queueRequest }) {
  const dispatch = useDispatch();

  const {
    submitted,
    isLoading,
    activeResultsTab,
    isError,
    gwas,
    LDFile,
    ldKey,
    select_qtls_samples,
    locus_quantification,
    locus_alignment,
    locus_table,
    isQueue,
  } = useSelector((state) => state.qtlsGWAS);

  const tabs = [
    {
      component: <LocusQC />,
      key: 'locus-qc',
      title: 'Locus QC',
      disabled: !submitted || isError,
    },
    {
      component: <LocusLD />,
      key: 'locus-ld',
      title: 'Locus LD',
      disabled:
        !submitted || isError || (!select_qtls_samples && !LDFile && !ldKey),
    },
    {
      component: <LocusAlignment />,
      key: 'locus-alignment',
      title: 'Locus Alignment',
      disabled: !submitted || isError,
    },
    {
      component: <LocusColocalization />,
      key: 'locus-colocalization',
      title: 'Locus Colocalization',
      disabled:
        !submitted ||
        isError ||
        !(gwas && gwas.data && Object.keys(gwas.data).length > 0),
    },
    {
      component: <LocusTable />,
      key: 'locus-table',
      title: 'Locus Table',
      disabled: !submitted || isError,
    },
    {
      component: <LocusQuantifiation />,
      key: 'locus-quantification',
      title: 'Locus Quantification',
      disabled:
        !submitted ||
        isError ||
        !(
          locus_quantification &&
          locus_quantification.data &&
          Object.keys(locus_quantification.data).length > 0
        ),
    },
    {
      component: <LocusDownload />,
      key: 'locus-download',
      title: 'Locus Download',
      disabled: !submitted || isError,
    },
  ];

  useEffect(() => {
    if (queueRequest) {
      dispatch(fetchResults({ request: queueRequest }));
    }
  }, [queueRequest]);

  return (
    <>
      <LoadingOverlay active={isLoading} />
      <Tabs
        // variant="pills"
        id="controlled-tab-example"
        activeKey={activeResultsTab}
        onSelect={(k) => {
          dispatch(updateQTLsGWAS({ activeResultsTab: k }));
        }}
      >
        {tabs.map((item, idx) => (
          <Tab
            key={item.key}
            eventKey={item.key}
            title={item.title}
            disabled={item.disabled}
            tabClassName={`mx-0 ${
              activeResultsTab == item.key ? 'font-weight-bold' : ''
            }`}
          >
            <div className="rounded-bottom border-left border-bottom border-right ">
              {(item.key == 'locus-alignment' || item.key == 'locus-table') &&
                (locus_alignment.data || locus_table.data) &&
                !isQueue &&
                submitted &&
                !isError && (
                  <>
                    <div className="px-3 py-2">
                      <QTLsGWASResultsForm />
                    </div>
                    <hr />
                  </>
                )}
              {item.component}
            </div>
          </Tab>
        ))}
      </Tabs>
    </>
  );
}
