    int getCubeX(long    m_volID)
    {
        if (m_volID >= 1e9)m_volID -= 1e9;
        return(((m_volID/10)%10000)/100);
    }

    int getCubeY(long    m_volID)
    {
        if (m_volID >= 1e9)m_volID -= 1e9;
        return((m_volID/10)%100);
    }

    int getCubeZ(long    m_volID)
    {
        if (m_volID >= 1e9)m_volID -= 1e9;
        return(m_volID/100000);
    }
